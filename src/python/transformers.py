import os
import sys
import subprocess
import time
import pickle
import shutil
import re
import numpy as np
import multiprocessing as mp
from datetime import timedelta
from setupparser import SetupParser
from copy import deepcopy
from gvar import VariableInteractionGraph
from gptlparser import gptl_parse_subset
from slicer import unslice
from reducer import reduce
from glob import glob

assert sys.version_info >= (3, 7), "Python version must be at least 3.7"

class ProseProjectTransformer:

    @staticmethod
    def load(path_to_transformer,resume=False):

        os.environ["PROSE_EXPERIMENT_DIR"] = os.getcwd()
        os.environ["PROSE_REPO_PATH"] = os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + "/..")

        if "Derecho" in os.path.basename(path_to_transformer):
            transformer = DerechoProseProjectTransformer._load(path_to_transformer)
        else:
            with open(path_to_transformer, "rb") as f:
                transformer = pickle.load(f)

        # if resuming, remove any directories of configurations that were not completely evaluated and reset
        if resume:
            for configuration_number in [int(dirname) for dirname in os.listdir(os.path.join(transformer.PROSE_EXPERIMENT_DIR, "prose_logs")) if dirname.isnumeric()]:
                if configuration_number > transformer.last_completed_configuration_number:
                    shutil.rmtree(os.path.join(transformer.PROSE_EXPERIMENT_DIR, "prose_logs/{:0>4}".format(configuration_number)))
            transformer.reset_project()

        return transformer


    def __new__(cls, *args, **kwargs):
        
        try:
            if cls != DerechoProseProjectTransformer and "Derecho" in SetupParser(args[0], working_dir=os.getcwd())._data:
                return DerechoProseProjectTransformer(args[0])
        except IndexError:
            pass

        return super(ProseProjectTransformer, cls).__new__(cls)

    def __init__(self, path_to_setup_file):
        print("** Setting Paths")

        os.environ["PROSE_REPO_PATH"] = os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + "/../..")
        os.environ["PROSE_EXPERIMENT_DIR"] = os.getcwd()
        self.PROSE_EXPERIMENT_DIR = os.environ["PROSE_EXPERIMENT_DIR"]

        self.SETUP = SetupParser(path_to_setup_file, working_dir=self.PROSE_EXPERIMENT_DIR)
        self.source_transformers_dict = {}
        self.G_proc = None
        self.GP_vertex_map = []
        self.G_var = None
        self.search_space = {}
        self.last_completed_configuration_number = -1
        self.timeout = -1

        try:
            self.ROSE_EXE_PATH = os.environ["ROSE_EXE_PATH"]
            self.PROSE_PLUGIN_PATH = os.environ["PROSE_PLUGIN_PATH"]
        except KeyError as e:
            raise Exception("[ERROR] ROSE_EXE_PATH or PROSE_PLUGIN_PATH environment variables not found.") from e

        print("** Creating required directories")
        try:
            if os.path.exists("prose_logs"):
                shutil.rmtree("prose_logs")
            if os.path.exists("prose_workspace"):
                shutil.rmtree("prose_workspace")
            os.makedirs("prose_logs")
            os.makedirs("prose_workspace")
            os.makedirs("prose_workspace/__search_progress")
            os.makedirs("prose_workspace/original_files")
            os.makedirs("prose_workspace/__profiling")
            os.makedirs("prose_workspace/__source_data")
            os.makedirs("prose_workspace/__timers")
            os.makedirs("prose_workspace/__timers/apply_configuration")
            os.makedirs("prose_workspace/__timers/compile")
            os.makedirs("prose_workspace/__timers/execute")
            os.makedirs("prose_workspace/__timers/evaluate")
            os.makedirs("prose_workspace/rmod_files")

        except:
            raise Exception("[ERROR] Failed to create required directories.")


    def _save(self):
        with open(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_workspace/__ProseProjectTransformer.pckl"), "wb") as f:
            pickle.dump(self, f)


    def preliminary_analysis(self):

        self._compile_plugin()

        # print("** Profiling application")
        # start = time.time()
        # self._conduct_profiling()
        # with open(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_logs/timers.txt"), "a") as f:
        #     f.write("{} profiling\n".format(timedelta(seconds=time.time()-start)))

        start = time.time()
        self._build_graphs()
        with open(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_logs/timers.txt"), "a") as f:
            f.write("{} building graphs\n".format(timedelta(seconds=time.time()-start)))

        start = time.time()
        self._construct_source_transformers()
        with open(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_logs/timers.txt"), "a") as f:
            f.write("{} constructing source transformers\n".format(timedelta(seconds=time.time()-start)))

        start = time.time()
        self._load_G_proc()
        self._load_G_var()
        with open(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_logs/timers.txt"), "a") as f:
            f.write("{} loading graphs\n".format(timedelta(seconds=time.time()-start)))

        start = time.time()
        self._propagate_constants()
        with open(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_logs/timers.txt"), "a") as f:
            f.write("{} propagating constants\n".format(timedelta(seconds=time.time()-start)))

        self._save()


    def generate_search_space(self, path_to_setup_file=None):

        if path_to_setup_file:
            os.environ["PROSE_EXPERIMENT_DIR"] = os.getcwd()
            self.PROSE_EXPERIMENT_DIR = os.environ["PROSE_EXPERIMENT_DIR"]
            self.SETUP = SetupParser(path_to_setup_file, working_dir=self.PROSE_EXPERIMENT_DIR)

        self._compile_plugin()

        start = time.time()
        self._reason_about_tuning_targets()
        with open(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_logs/timers.txt"), "a") as f:
            f.write("{} generating search space\n".format(timedelta(seconds=time.time()-start)))

        print("** compiling and caching for future incremental builds")
        
        subprocess.run(
            self.SETUP['build']['cmd'],
            check=True,
            env=os.environ.copy(),
            cwd=self.SETUP['build']['working_dir'],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            executable="/bin/bash",
            shell=True,
        )

        return deepcopy(self.search_space)


    def search(self, search_algorithm):
        start = time.time()

        if float(self.SETUP['eval']['cost_threshold']) > 0:
            search_algorithm.set_cost_threshold(float(self.SETUP['eval']['cost_threshold']))

        configuration_number = search_algorithm.completed_config_counter
        configuration_dict = search_algorithm.get_next()
        while configuration_dict:
            total_cost, targeted_subset_cost = self._test_configuration(configuration_dict=configuration_dict, configuration_number=configuration_number)
            if targeted_subset_cost != 0:
                cost = targeted_subset_cost
            else:
                cost = total_cost

            configuration_dict["cost"] = cost
            if self.timeout < 0:
                if self.SETUP['run']['timeout'] == '0':
                    self.timeout = int(np.ceil(total_cost * 3.0))
                else:
                    self.timeout = int(float(self.SETUP['run']['timeout']))

            search_algorithm.feedback(configuration_dict)
            self.last_completed_configuration_number = configuration_number            
            self._save()
            
            configuration_dict = search_algorithm.get_next()
            configuration_number += 1

        with open(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_logs/timers.txt"), "a") as f:
            f.write("{} search\n".format(timedelta(seconds=time.time()-start)))

        self._save()
        self.report(final=True)


    def test_configuration(self, configuration_file_path, path_to_setup_file=None, node_name=""):

        if path_to_setup_file:
            os.environ["PROSE_EXPERIMENT_DIR"] = os.getcwd()
            self.PROSE_EXPERIMENT_DIR = os.environ["PROSE_EXPERIMENT_DIR"]
            self.SETUP = SetupParser(path_to_setup_file, working_dir=self.PROSE_EXPERIMENT_DIR)

        assert(os.path.exists(configuration_file_path))
        assert(os.getcwd() == self.PROSE_EXPERIMENT_DIR)
        self.reset_project()
        try:
            shutil.copyfile(configuration_file_path, os.path.join(os.path.dirname(configuration_file_path), "config"))
        except shutil.SameFileError:
            pass

        print("\n============ Testing Configuration ============\n")

        working_dir = self.PROSE_EXPERIMENT_DIR
        configuration_dir = os.path.dirname(configuration_file_path)

        self._compile_plugin()
        flag = self._apply_configuration(working_dir, configuration_dir)
        if flag < 0:
            print("\t !! unable to apply transformations to source code")

        else:

            flag = self._compile(working_dir, configuration_dir)
            if flag < 0:
                print("\t !! unable to compile transformed source code")
            else:

                exception = self._execute(working_dir, configuration_dir, node_name)
                if exception:
                    raise exception
                else:

                    total_cost, targeted_subset_cost = self._evaluate(working_dir, configuration_dir)
                    if targeted_subset_cost != 0:
                        cost = targeted_subset_cost
                    else:
                        cost = total_cost

                    if np.isinf(targeted_subset_cost):
                        print("\t **[INVALID] gptl parsing error")
                    elif cost > 0:
                        if targeted_subset_cost != 0:
                            print("\t ** [PASSED] subset cost = {:.3f} (total cost = {:.3f})".format(targeted_subset_cost, total_cost))
                        else:
                            print("\t ** [PASSED] cost = {:.3f}".format(cost))
                    elif cost < 0:
                        if targeted_subset_cost != 0:
                            print("\t ** [FAILED] subset cost = {:.3f} (total cost = {:.3f}) but error threshold was exceeded".format(np.abs(targeted_subset_cost), np.abs(total_cost)))
                        else:
                            print("\t ** [FAILED] cost = {:.3f} but error threshold was exceeded".format(np.abs(cost)))

        os.remove(os.path.join(os.path.dirname(configuration_file_path), "config"))
        self.reset_project()


    def _test_configuration(self, configuration_dict, configuration_number, working_dir=".", node_name=""):

        print("\n============ Testing Configuration {} ============\n".format(configuration_number))
        
        configuration_dir = os.path.join(working_dir, "prose_logs/{:0>4}".format(configuration_number))
        configuration_number = "{}".format(configuration_number)
        os.makedirs(configuration_dir)
        
        # write out the configuration to be read in by the ROSE plugin
        with open(os.path.join(configuration_dir, "config"), "w") as f:
            for var_name, varKind in configuration_dict["config"].items():
                f.write("{},{}\n".format(var_name, varKind))

        total_cost = -1
        targeted_subset_cost = -1
        
        if self.SETUP['run']['execution_filtering'].lower() == "true":
            predicted_cost_ratio = self.G_var.get_cost_ratio(configuration_dict)
        else:
            predicted_cost_ratio = 1

        if predicted_cost_ratio < 1: 
            message = "{:0>4}: [FAILED] (cost model) predicted cost ratio = {}".format(configuration_number, predicted_cost_ratio)
            os.rename(os.path.join(configuration_dir, "config"), os.path.join(configuration_dir, "config_FAILED"))
            total_cost = np.inf
            targeted_subset_cost = np.inf
        else:

            start = time.time()
            flag = self._apply_configuration(working_dir, configuration_dir, node_name)
            os.mknod(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_workspace/__timers/apply_configuration/{:0>4}_{}".format(configuration_number, time.time()-start)))

            if flag < 0:
                message = "{:0>4}: [INVALID] (plugin error) unable to apply transformations to source code".format(configuration_number)
                os.rename(os.path.join(configuration_dir, "config"), os.path.join(configuration_dir, "config_INVALID"))
                total_cost = np.inf
                targeted_subset_cost = np.inf
            else:

                start = time.time()
                flag = self._compile(working_dir, configuration_dir, node_name)
                os.mknod(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_workspace/__timers/compile/{:0>4}_{}".format(configuration_number, time.time()-start)))

                if flag < 0:
                    message = "{:0>4}: [INVALID] (compilation error) unable to compile transformed source code".format(configuration_number)
                    os.rename(os.path.join(configuration_dir, "config"), os.path.join(configuration_dir, "config_INVALID"))
                    total_cost = np.inf
                    targeted_subset_cost = np.inf
                else:

                    start = time.time()
                    exception = self._execute(working_dir, configuration_dir, node_name)
                    os.mknod(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_workspace/__timers/execute/{:0>4}_{}".format(configuration_number, time.time()-start)))

                    if exception:
                        total_cost = np.inf
                        targeted_subset_cost = np.inf
                        if isinstance(exception, subprocess.TimeoutExpired):
                            message = "{:0>4}: [FAILED] (timeout) timeout {} exceeded".format(configuration_number, self.timeout)
                            os.rename(os.path.join(configuration_dir, "config"), os.path.join(configuration_dir, "config_FAILED"))
                        else:
                            message = "{:0>4}: [FAILED] (runtime failure) unable to execute transformed source code (see outlog.txt for details)".format(configuration_number)
                            os.rename(os.path.join(configuration_dir, "config"), os.path.join(configuration_dir, "config_FAILED"))
                    else:

                        start = time.time()  
                        total_cost, targeted_subset_cost = self._evaluate(working_dir, configuration_dir)

                        os.mknod(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_workspace/__timers/evaluate/{:0>4}_{}".format(configuration_number, time.time()-start)))

                        # if _evaluate returned 0, it is because it couldn't cast the return from the eval command
                        # as a float which we assume means that the program failed gracefully but failed nonetheless, e.g., returning NaN
                        if total_cost == 0:
                            message = "{:0>4}: [FAILED] (runtime failure) unable to execute transformed source code (see outlog.txt for details)".format(configuration_number)
                            os.rename(os.path.join(configuration_dir, "config"), os.path.join(configuration_dir, "config_FAILED"))
                            total_cost = np.inf
                            targeted_subset_cost = np.inf
                        elif np.isinf(targeted_subset_cost):
                            message = "{:0>4}: [INVALID] gptl parsing error".format(configuration_number)
                            os.rename(os.path.join(configuration_dir, "config"), os.path.join(configuration_dir, "config_INVALID"))
                            total_cost = np.inf
                            targeted_subset_cost = np.inf
                        
                        # otherwise, the program terminated gracefully and has an associated cost as determined by the provided eval program
                        else:

                            if targeted_subset_cost != 0:
                                cost = targeted_subset_cost
                            else:
                                cost = total_cost

                            os.mknod(os.path.join(configuration_dir, f"COST_{np.abs(cost)}"))
    
                            if cost > 0:
                                if targeted_subset_cost != 0:
                                    message = "{:0>4}: [PASSED] subset cost = {:.3f} (total cost = {:.3f})".format(configuration_number, targeted_subset_cost, total_cost)
                                else:
                                    message = "{:0>4}: [PASSED] cost = {:.3f}".format(configuration_number, cost)
                                os.rename(os.path.join(configuration_dir, "config"), os.path.join(configuration_dir, "config_PASSED"))
                            elif cost < 0:
                                if targeted_subset_cost != 0:
                                    message = "{:0>4}: [FAILED] subset cost = {:.3f} (total cost = {:.3f}) but error threshold was exceeded".format(configuration_number, np.abs(targeted_subset_cost), np.abs(total_cost))
                                else:
                                    message = "{:0>4}: [FAILED] cost = {:.3f} but error threshold was exceeded".format(configuration_number, np.abs(cost))
                                os.rename(os.path.join(configuration_dir, "config"), os.path.join(configuration_dir, "config_FAILED"))
                                total_cost = np.inf
                                targeted_subset_cost = np.inf

        # there has to be something wrong with our transformer if the first config is already invalid.
        if configuration_number == 0:
            assert(total_cost > 0)

        print("\t **{}".format(message[message.find(":") + 1:]))
        with open(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_workspace/__search_progress/{:0>4}".format(configuration_number)), "w") as f:
            f.write("{}\n".format(message))

        self.reset_project()

        return total_cost, targeted_subset_cost


    def report(self, final=False):

        # search log processing
        buffer = []
        failed_costs = []
        passed_costs = []
        for message_file in sorted(os.listdir(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_workspace/__search_progress")), key=lambda x: int(x)):
            with open(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_workspace/__search_progress", message_file), "r") as f:
                
                line = f.readline()

                tokens = line.split()
                if "cost =" in line and "[PASSED]" in tokens:
                    cost = float(tokens[tokens.index("=") + 1])
                    failed_costs.append(np.inf)
                    passed_costs.append(cost)
                elif "cost =" in line and "[FAILED]" in tokens:
                    cost = float(tokens[tokens.index("=") + 1])
                    failed_costs.append(cost)
                    passed_costs.append(np.inf)
                else:
                    failed_costs.append(np.inf)
                    passed_costs.append(np.inf)
                    
                buffer.append(line)

        buffer.append("\n          Original cost = {:.3f}\n".format(passed_costs[0]))
        buffer.append("Best FAILED: {:0>4}, cost = {:.3f} ({:.3f}x speedup)\n".format(np.argmin(failed_costs), np.min(failed_costs), passed_costs[0]/np.min(failed_costs)))
        buffer.append("Best PASSED: {:0>4}, cost = {:.3f} ({:.3f}x speedup)\n".format(np.argmin(passed_costs), np.min(passed_costs), passed_costs[0]/np.min(passed_costs)))
        
        if not final:
            for line in buffer:
                print(line, end="")
            return
        else:
            with open(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_logs/search_log.txt"), "w") as f:
                f.writelines(buffer)

            # timer processing
            for timer_label in ["apply_configuration", "compile", "execute", "evaluate"]:
                buffer = []
                for message_file in sorted(os.listdir(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_workspace/__timers/{}".format(timer_label)))):
                    _, time = message_file.split("_")
                    buffer.append(float(time))
                with open(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_logs/timers.txt"), "a") as f:
                    f.write("\n{} details:\n".format(timer_label))
                    f.write("\t    total: {}\n".format(timedelta(seconds=np.sum(buffer))))
                    f.write("\t     mean: {}\n".format(timedelta(seconds=np.mean(buffer))))
                    f.write("\t variance: {}\n".format(timedelta(seconds=np.var(buffer))))


    def reset_project(self):
        for _, transformer in self.source_transformers_dict.items():
            transformer.reset()


    def _compile_plugin(self):
        # make sure plugin is compiled
        try:
            subprocess.run(
                "make",
                check=True,
                env=os.environ.copy(),
                cwd=self.PROSE_PLUGIN_PATH,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                executable="/bin/bash",
                shell=True,
            )
        except:
            raise Exception("[ERROR] Plugin failed to compile.")

    def _build_graphs(self):

        os.makedirs("prose_workspace/temp_graph")

        self.src_to_transform = reduce(self.SETUP)

        # link graphs into one single graph
        prose_command = [
            "{}/rose-compiler".format(self.ROSE_EXE_PATH),
            "-rose:plugin_lib","{}/ProsePlugin.so".format(self.PROSE_PLUGIN_PATH),
            "-rose:plugin_action","prose-link-graph",
            "-rose:plugin_arg_prose-link-graph", "prose_workspace/",
            "-rose:skip_syntax_check",
            "-rose:skip_unparse",
            "-rose:skipfinalCompileStep",
            "-rose:output", "prose_workspace/__temp",
            "-Drose_comp",
            "-I{}/prose_workspace/rmod_files".format(self.PROSE_EXPERIMENT_DIR),
            self.SETUP['target']['additional_plugin_flags'],
            os.path.abspath(self.SETUP['target']['src_files'].split("|")[0])
        ]
        subprocess.run(
            ' '.join(prose_command),
            check=True,
            env=os.environ.copy(),
            executable="/bin/bash",
            shell=True,
        )

        # clean intermediate files
        shutil.rmtree("prose_workspace/temp_graph")
        os.remove("prose_workspace/__temp")


    def _load_G_proc(self):

        # Parse G_proc.dot
        with open("prose_logs/__G_proc.dot", "r") as f:
            lines = f.readlines()

            # load all vertices, already in vertexD order
            for line in lines:
                if "{" in line or "}" in line:
                    continue    # skip first and last line of .dot file
                elif line.find("->") == -1:
                    scope_name = line.split('"')[1]
                    assert( scope_name == scope_name.lower() )
                    self.GP_vertex_map.append(scope_name)

            # initialize G_proc
            self.G_proc = np.zeros((len(self.GP_vertex_map), len(self.GP_vertex_map)))

            # load all edges
            for line in lines:
                if line.find("->") != -1:
                    source_vertexD = int(line.split('-')[0].strip())
                    target_vertexD = int(line.split('[')[0].split('>')[1].strip())
                    self.G_proc[source_vertexD][target_vertexD] += 1


    def _construct_source_transformers(self):
        # construct source transformers for all source files in
        # project and load all scope_name to source transformer
        # mappings; the format of the datafiles is: first line is the
        # path to the source in the project second line is all of the
        # scopes in this source separated by + signs

        for path_to_data_file in [os.path.join("prose_workspace/__source_data/", x) for x in os.listdir("prose_workspace/__source_data/")]:
            transformer = ProseSourceTransformer(path_to_data_file, experiment_dir=self.PROSE_EXPERIMENT_DIR)
            self.source_transformers_dict[transformer.get_name()] = transformer

            for scope_name in transformer.get_scope_names():
                self.source_transformers_dict[scope_name] = transformer

        # reset project after graph building (i.e., reset sliced files to original src files)
        self.reset_project()


    def _apply_configuration(self, working_dir, configuration_dir, node_name=""):

        # write out the list of files to be transformed; will be read in by the plugin
        with open(os.path.join(working_dir, "prose_workspace/__target_files.txt"), "w") as f:
            for src_file_path in self.src_to_transform:
                f.write(src_file_path + "\n")
                self.source_transformers_dict[os.path.basename(src_file_path[:src_file_path.rfind(".")]).lower()].pre_transform_process(working_dir)

        prose_command = [
            "{}/rose-compiler".format(self.ROSE_EXE_PATH),
            "-rose:plugin_lib","{}/ProsePlugin.so".format(self.PROSE_PLUGIN_PATH),
            "-rose:plugin_action","prose-apply-configuration",
            "-rose:plugin_arg_prose-apply-configuration", working_dir,
            "-rose:plugin_arg_prose-apply-configuration", configuration_dir,
            "-rose:skip_syntax_check",
            "-rose:skip_unparse",
            "-rose:skipfinalCompileStep",
            "-Drose_comp",
            "-I{}/prose_workspace/rmod_files".format(self.PROSE_EXPERIMENT_DIR),
            self.SETUP['target']['additional_plugin_flags'],
            os.path.join(working_dir, self.SETUP['target']['src_files'].split("|")[0])
        ]

        if node_name:
            prose_command = [f"ssh {node_name} 'source {self.SETUP['Derecho']['env_script']} && cd {working_dir} &&"] + prose_command + ["'"]
        prose_command = ' '.join(prose_command)

        try:
            subprocess.run(
                prose_command,
                check = True,
                stdout = subprocess.DEVNULL,
                stderr = subprocess.DEVNULL,
                cwd = working_dir,
                env = os.environ.copy(),
                executable="/bin/bash",
                shell=True,
            )
        except subprocess.CalledProcessError as e:
            return -1
   
        return 0


    def _compile(self, working_dir, configuration_dir, node_name=""):

        # unslice and move all of the transformed source code to the proper
        # locations in order to be compiled
        print("\t ** unslicing")
        with open(os.path.join(working_dir, "prose_workspace/__target_files.txt"), "r") as f:
            target_file_paths = f.readlines()
            for target_file_path in target_file_paths:
                transformer = self.source_transformers_dict[os.path.basename(target_file_path[:target_file_path.rfind(".")]).lower()]
                transformer.post_transform_process(working_dir, configuration_dir, self.SETUP)

        command = self.SETUP['build']['partial_build_cmd'].split()

        # used for one-off configurations on Derecho being tested outside of the context of a full precision tuning run
        # see scripts/prose_test_single_configuration.py
        if node_name == "offline":
            command = [f"source {self.SETUP['Derecho']['env_script']} && "] + command

        elif node_name:
            command = [f"ssh {node_name} 'source {self.SETUP['Derecho']['env_script']} && cd {os.path.join(working_dir, self.SETUP['build']['working_dir'])} &&"] + command + ["'"]
        command = " ".join(command)

        print("\t ** compiling")
        try:
            subprocess.run(
                command,
                check = True,
                shell = True,
                executable="/bin/bash",
                env = os.environ.copy(),
                cwd = os.path.join(working_dir, self.SETUP['build']['working_dir']),
                stdout = subprocess.PIPE,
                stderr = subprocess.STDOUT,
            )

        except subprocess.CalledProcessError as e:
            compile_log = e.output.decode("utf-8").splitlines()
            with open(os.path.join(configuration_dir, "compile_error.txt"), "a") as f:
                for line in compile_log:
                    f.write(line+"\n")
            return -2

        return 0


    def _execute(self, working_dir, configuration_dir, node_name=""):

        print("\t ** executing")

        if self.timeout > 0:
            timeout = self.timeout
        else:
            timeout = None

        command = self.SETUP['run']['cmd'].split()
        if node_name:
            
            # used for one-off configurations on Derecho being tested outside of the context of a full precision tuning run
            # see scripts/prose_test_single_configuration.py
            if node_name == "offline":
                command = [f"qcmd -- 'source {self.SETUP['Derecho']['env_script']} && "] + command + ["'"]
                timeout = None

            # otherwise, the node name is one of those reserved via PBS
            else:

                # insert the name of the node to be used under this mpi invocation
                for mpi_keyword in ["mpiexec", "mpirun"]:
                    try:
                        command = command[:command.index(mpi_keyword) + 1] + ["--hosts", node_name] + command[command.index(mpi_keyword) + 1:]
                        break
                    except ValueError:
                        pass

                # if we didn't successfully add the node_name (because the run command does not have mpiexec or mpirun in it), append the node name
                if node_name not in command:
                    command = command + [node_name]

                # compose the command to be executed over ssh on the specified node
                # source the env script, cd to the correct directory, and execute
                command = [f"ssh {node_name} 'source {self.SETUP['Derecho']['env_script']} && cd {os.path.join(working_dir, self.SETUP['run']['working_dir'])} &&"] + command + ["'"]

        command = " ".join(command)

        try:
            with open(os.path.join(configuration_dir, "outlog.txt"), "w") as outfile:
                subprocess.run(
                            command,
                            check=True,
                            env=os.environ.copy(),
                            cwd=os.path.join(working_dir, self.SETUP['run']['working_dir']),
                            stdout=outfile,
                            stderr=outfile,
                            timeout=timeout,
                            shell=True,
                            executable="/bin/bash"
                        )      
            
            return None

        except subprocess.SubprocessError as e:
            return e
        

    def _evaluate(self, working_dir, configuration_dir):

        print("\t ** evaluating")

        command = self.SETUP['eval']['cmd'].split()
        if self.SETUP['eval']['pass_log_path'].strip().lower() == "true":
            original_configuration_dir = os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_logs/0000")
            command += [configuration_dir, original_configuration_dir]            
        command = " ".join(command)

        evaluation_output = subprocess.check_output(
                        command,
                        env=os.environ.copy(),
                        cwd=os.path.join(working_dir, self.SETUP['eval']['working_dir']),
                        shell=True,
                        executable="/bin/bash",
                    )
        
        evaluation_output = evaluation_output.decode("utf-8").splitlines()
        try:
            total_cost = float(evaluation_output[0])
        except ValueError:
            return 0, 0

        if len(glob(os.path.join(working_dir, "timing.*"))) == 0:
            targeted_subset_cost = 0.0
        else:

            # parse gptl subset timing if available
            if self.SETUP['target']['ignore_patterns']:
                targeted_subset_cost = gptl_parse_subset(self.SETUP['target']['search_patterns'].split("|"), self.SETUP['target']['ignore_patterns'].split("|"), working_dir)
            else:
                targeted_subset_cost = gptl_parse_subset(self.SETUP['target']['search_patterns'].split("|"), [], working_dir)
            
            # save the info generated from the call above
            subprocess.run(f"tar -czf {configuration_dir}/gptl_timing.tar.gz timing.* && rm timing.*", shell=True, executable="/bin/bash", cwd=working_dir)
            os.rename(os.path.join(working_dir, "gptl_subset_info.pckl"), os.path.join(configuration_dir, "gptl_subset_info.pckl"))
            
            # if the total cost is negative, this means it failed the error check
            # flip the sign of the subset cost as well
            if total_cost < 0:
                targeted_subset_cost = -1 * targeted_subset_cost

                # signal mismatching signs via inf; this means that gptl parsing failed and returned a negative cost
                if total_cost * targeted_subset_cost < 0:
                    targeted_subset_cost = np.inf

        return total_cost, targeted_subset_cost
    

    def _is_target(self, scoped_name):

        targeted = 0

        for search_pattern in self.SETUP['target']['search_patterns'].split("|"):
            if search_pattern and re.search(r"{}".format(search_pattern), scoped_name, re.IGNORECASE):
                targeted = 1
                break

        for ignore_pattern in self.SETUP['target']['ignore_patterns'].split("|"):
            if ignore_pattern and re.match(r"{}".format(ignore_pattern), scoped_name, re.IGNORECASE):
                targeted = -1
                break

        return targeted


    def _reason_about_tuning_targets(self):

        possible_kinds = [4,8,16]
        targeted_vars = []
        for var_name in self.G_var.var_names:
            if self._is_target(var_name) > 0: # matches a search pattern, doesn't match any ignore patterns
                targeted_vars.append(var_name)

        if self._is_target("::MOM_continuity_PPM::foo") > 0:
            targeted_vars = sorted(targeted_vars, key=lambda x : x.split(":")[4].replace("_y", "_y1").replace("_x","_x1").replace("meridional", "merid").replace("zonal", "").replace("merid","")[::-1])

        for var_name in targeted_vars:
            scope_name = var_name[:var_name.rfind("::")]
            self.search_space[var_name] = [str(x) for x in possible_kinds if x <= self.source_transformers_dict[scope_name].variable_profile[var_name]['kind']] 

        # write out ignored scopes to be read in by plugin
        with open("prose_workspace/ignore_scopes.txt", "w") as f:
            for scope_name in self.GP_vertex_map:
                if self._is_target(scope_name) < 0:
                    f.write(scope_name + "\n")


    def _load_G_var(self):
        self.G_var = VariableInteractionGraph("prose_logs/__G_var.dot", self.source_transformers_dict)


    def _propagate_constants(self):

        def flows_to(var1, var2):
            scope_name1 = var1[:var1.rfind("::")]
            scope_name2 = var2[:var2.rfind("::")]

            failed = False
            try:
                scope1_vertexD = self.GP_vertex_map.index(scope_name1)
            except:
                print("Couldn't find scope for {}; likely an imported constant? skipping".format(var1))
                failed = True
            try:
                scope2_vertexD = self.GP_vertex_map.index(scope_name2)
            except:
                print("Couldn't find scope for {}; likely an imported constant? skipping".format(var2))
                failed = True
            if failed:
                return False

            if self.G_proc[scope1_vertexD][scope2_vertexD] > 0:
                return True
            else:
                return False

        inter_bindings = {}
        # build up a directed graph (adjacency list) of
        # interprocedural variable bindings
        with open("prose_workspace/__inter_bound_variables.txt", 'r') as f:
            for line in f:
                binding = [x for x in line.strip().split(";") if x != ""]
                for var_name in binding:
                    inter_bindings[var_name] = inter_bindings.get(var_name, set()).union(set([x for x in binding if flows_to(var_name, x)]))

        discovered_constant_vars = set()
        constant_bindings = []

        # read constant list file, each line is one const
        with open("prose_workspace/constant_list.txt", 'r') as f:
            for line in f:
                var_name = line
                discovered_constant_vars.add(var_name.strip())

        # propagate constants
        for constant_var in discovered_constant_vars:
            stack = [constant_var]
            visited = set()
            const_propagation = set()

            while len(stack) != 0:
                from_var = stack.pop()
                if from_var not in visited and from_var in inter_bindings:
                    visited.add(from_var)
                    const_propagation.add(from_var)
                    stack = stack + list(inter_bindings[from_var])

            const_propagation.discard(constant_var)
            constant_bindings.append(const_propagation)

        # write out propagated constant list
        with open("prose_workspace/constant_list.txt", 'w') as f:
            for binding in constant_bindings:
                for var_name in binding:
                    f.write(var_name + "\n")


class DerechoProseProjectTransformer(ProseProjectTransformer):

    @staticmethod
    def _load(path_to_transformer):
        with open(path_to_transformer, "rb") as f:
            return pickle.load(f)

    def __new__(cls, *args, **kwargs):
        return super(DerechoProseProjectTransformer, cls).__new__(cls)

    def __init__(self, path_to_setup_file):

        super(DerechoProseProjectTransformer, self).__init__(path_to_setup_file)
        
        self.scratch_path = os.path.join("/glade/derecho/scratch", os.environ['USER']) + self.SETUP['machine']['project_root']
        if os.path.exists(self.scratch_path):
            shutil.rmtree(self.scratch_path)
            try:
                os.removedirs(os.path.abspath(os.path.join(self.scratch_path, os.pardir)))
            except OSError:
                pass
        os.makedirs(self.scratch_path)


    def _save(self):
        with open(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_workspace/__DerechoProseProjectTransformer.pckl"), "wb") as f:
            pickle.dump(self, f)
    

    def generate_search_space(self, path_to_setup_file=None):

        search_space = super().generate_search_space(path_to_setup_file)
    
        if os.path.exists(os.path.join(self.scratch_path, "original")):
            shutil.rmtree(os.path.join(self.scratch_path, "original"))

        shutil.copytree(
            self.SETUP['machine']['project_root'],
            os.path.join(self.scratch_path, "original"),
            ignore=shutil.ignore_patterns(*(self.SETUP['Derecho']['copy_ignore'].split("|") + ["prose_logs", "prose_workspace", ".git"])),
            symlinks=True,
        )
        shutil.copytree(
            os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_workspace"),
            os.path.join(self.scratch_path, "original", os.path.relpath(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_workspace"), start=self.SETUP['machine']['project_root'])),
            ignore=shutil.ignore_patterns("__*"),
            symlinks=True,
        )
        shutil.copytree(
            os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_logs"),
            os.path.join(self.scratch_path, "original", os.path.relpath(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_logs"), start=self.SETUP['machine']['project_root'])),
            ignore=shutil.ignore_patterns("__*"),
            symlinks=True,
        )

        return search_space


    def search(self, search_algorithm):
        start = time.time()

        if float(self.SETUP['eval']['cost_threshold']) > 0:
            search_algorithm.set_cost_threshold(float(self.SETUP['eval']['cost_threshold']))

        result = subprocess.run(
            'for x in $(cat $PBS_NODEFILE | sort | uniq | cut -d"". -f1); do if [[ "$(hostname)" == ${x} ]]; then continue; else echo ${x}; fi; done',
            shell=True,
            stdout=subprocess.PIPE,
            text=True,
            executable="/bin/bash"
        )
        node_list = [x for x in result.stdout.split("\n") if x]
        shared_q_avail_nodes = mp.Queue(maxsize=len(node_list))
        for node_name in node_list:
            shared_q_avail_nodes.put(node_name)

        batch = search_algorithm.get_next_batch()
        while batch:

            shared_q_process_returns = mp.Queue(maxsize=len(batch))
            for i in range(len(batch)):
                node_name = shared_q_avail_nodes.get(block=True)
                mp.Process(
                    target=self._test_configuration,
                    args=(
                        batch[i],
                        self.last_completed_configuration_number + 1 + i,
                        node_name,
                        shared_q_avail_nodes,
                        shared_q_process_returns,
                    )
                ).start()

            costs = [(None, None)]*len(batch)
            for i in range(len(batch)):
                configuration_number, cost = shared_q_process_returns.get(block=True)
                costs[configuration_number - self.last_completed_configuration_number - 1] = cost

            for i in range(len(batch)):
                if self.timeout < 0: # should only happen on the default (first) configuration
                    if self.SETUP['run']['timeout'] == '0':
                        try:
                            self.timeout = int(np.ceil(costs[0][0] * 3.0))
                        except OverflowError:
                            print("Default configuration appears to have failed; exiting.")
                            exit(1)
                    else:
                        self.timeout = int(float(self.SETUP['run']['timeout']))

                # take subset cost if present; otherwise, take total cost
                if costs[i][1] != 0:
                    batch[i]["cost"] = costs[i][1]
                else:
                    batch[i]["cost"] = costs[i][0]

                search_algorithm.feedback(batch[i])
            
            self.last_completed_configuration_number += len(batch)
            batch = search_algorithm.get_next_batch()

            self._save()

        with open(os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_logs/timers.txt"), "a") as f:
            f.write("{} search\n".format(timedelta(seconds=time.time()-start)))

        shutil.rmtree(self.scratch_path, ignore_errors=True)
        try:
            os.removedirs(os.path.abspath(os.path.join(self.scratch_path, os.pardir)))
        except OSError:
            pass
        
        self._save()
        self.report(final=True)

    def test_configuration(self, configuration_file_path, path_to_setup_file=None):
        super().test_configuration(configuration_file_path, path_to_setup_file, node_name="offline")


    def _test_configuration(self, configuration_dict, configuration_number, node_name, shared_q_avail_nodes, shared_q_process_returns):

        # copy backed-up "original" directory in the scratch space into its own experiment directory
        new_project_root = os.path.join(self.scratch_path, "{:0>4}".format(configuration_number))
        if os.path.exists(new_project_root):
            shutil.rmtree(new_project_root)
        shutil.copytree(
            os.path.join(self.scratch_path, "original"),
            new_project_root,
            symlinks=True,
        )
        
        working_dir = os.path.abspath(os.path.join(new_project_root, os.path.relpath(self.PROSE_EXPERIMENT_DIR, start=self.SETUP['machine']['project_root'])))

        cost = super()._test_configuration(configuration_dict, configuration_number, working_dir, node_name)
        shared_q_avail_nodes.put(node_name)
        shared_q_process_returns.put((configuration_number, cost))

        # copy everything back to the OG experiment directory and clean up 
        shutil.copytree(
            os.path.join(new_project_root, os.path.join(os.path.relpath("prose_logs", start=self.SETUP['machine']['project_root']),"{:0>4}".format(configuration_number))),
            os.path.join(self.PROSE_EXPERIMENT_DIR, "prose_logs/{:0>4}".format(configuration_number))
        )
        shutil.rmtree(new_project_root)


class ProseSourceTransformer:

    def __init__(self, path_to_data_file, experiment_dir):

        self.name = os.path.basename(path_to_data_file).lower()
        with open(path_to_data_file, "r") as f:
            lines = f.readlines()

        self.path_original_file = os.path.relpath(lines[0].strip(), start=experiment_dir)
        self.name_original_file = os.path.basename(self.path_original_file)
        self.scope_names = [ x for x in lines[1].strip().split("+") if x != "" ]
        for scope_name in self.scope_names:
            assert( scope_name == scope_name.lower() )

        self.variable_profile = {}

        for line in lines[2:]:
            var_name, var_type, var_kind = tuple(line.split(","))

            assert( var_name == var_name.lower() )

            var_type = var_type.split("=")[-1]
            if "(" in var_type:
                var_type, var_dim = tuple(var_type.split("("))
                var_dim = float(var_dim[:-1])
            else:
                var_dim = 1
            var_kind = int(var_kind.split("=")[-1])
            self.variable_profile[var_name] = {
                "type" : var_type,
                "dim"  : var_dim,
                "kind" : var_kind,
            }

        self.reset()


    def get_name(self):
        return self.name


    def get_scope_names(self):
        return self.scope_names


    def pre_transform_process(self, working_dir):
        shutil.copy("prose_workspace/original_files/{}.slice".format(self.name_original_file), os.path.join(working_dir, self.path_original_file))

    def post_transform_process(self, working_dir, configuration_dir, SETUP):

        # case insensitive search for transformed file
        matches = [ x for x in os.listdir(configuration_dir) if x.lower() == self.name_original_file.lower() ]
        assert(len(matches) == 1)

        # rename to be sure that file name case matches what is expected
        transformed_file_name = matches[0]
        os.rename(os.path.join(configuration_dir, transformed_file_name), os.path.join(configuration_dir, self.name_original_file))

        unslice(os.path.join(configuration_dir, self.name_original_file), "prose_workspace/original_files/{}.orig".format(self.name_original_file), SETUP)

        # move transformed file to its original location
        shutil.copy(
            os.path.join(configuration_dir, self.name_original_file),
            os.path.join(working_dir, self.path_original_file)
        )

    def reset(self):
        shutil.copy("prose_workspace/original_files/{}.orig".format(self.name_original_file), self.path_original_file)