#!/usr/bin/env python3
import os
import re
from glob import glob
import numpy as np
import sys
from copy import deepcopy
import pickle

class GPTLParseFail(Exception):
    pass

def wrapper_procedure(name):
    if re.search(r"[a-z0-9_](_wrap_[0-9][0-9][0-9]|_wrapper(_id[0-9][0-9][0-9])?_[0-9a-gx]+_to_[0-9a-gx])", name, re.IGNORECASE):
        return True
    else:
        return False

def gptl_parse_subset(search_patterns, ignore_patterns, working_dir=".", pprint=False):

    targeted_subset_runtimes = {}
    for gptl_timing_file in glob(os.path.join(working_dir, "timing.*")):
        
        try:
            with open(gptl_timing_file, "r", errors='replace') as f:
                lines = [line.replace(") [0x","__[0x", 1) for line in f.readlines()]

            processed = set()

            # find idx for start of subset timing information
            subset_timing_idx = 0
            while subset_timing_idx + 1 < len(lines) and not lines[subset_timing_idx].strip().startswith("Stats for thread 0:"):
                subset_timing_idx += 1
            subset_timing_idx += 1

            # pass over subset timing information to gather the index from which to start gathering runtimes in each line
            # and taking note of multi-parent entries and their runtimes
            # also, check if the gptl log file mis-printed; if so, we will skip this particular log
            i = subset_timing_idx
            line_length_invariant = len(lines[subset_timing_idx])
            while i + 1 < len(lines) and not lines[i + 1].startswith("Overhead sum ="):
                i += 1
                if len(lines[i]) != line_length_invariant:
                    raise GPTLParseFail(f"Parsing of {gptl_timing_file} failed on the following line:\n{lines[i]}")

            # gather long-name-translation information
            long_name_to_alias_map = {}
            alias_to_long_name_map = {}
            while i + 1 < len(lines) and not lines[i].strip().startswith("thread 0 long name translations"):
                i += 1
            while i + 1 < len(lines) and not lines[i + 1].strip().startswith("Multiple parent info for thread 0:"):
                i += 1
                if " = " in lines[i]: 
                    alias, long_name = (x.strip() for x in lines[i].split("=", maxsplit=1))
                    long_name_to_alias_map[long_name] = alias
                    alias_to_long_name_map[alias] = long_name

            # gather target subset runtime if it exists
            t = {}
            i = subset_timing_idx
            while i + 1 < len(lines) and not lines[i+1].startswith("Overhead sum ="):
                i += 1

                time = float(lines[i][2:].split()[3])
                name = lines[i][2:].split()[0]
                if name in alias_to_long_name_map:
                    name = alias_to_long_name_map[name]
                name = name.lower()

                # targeted procedures will have :: in their timer name;
                # ones that don't are wrappers
                if "::" not in name:

                    # special exception for mpas wrapper routines for wrappers around targeted work routines in atm_time_integration;
                    # manual inspection says that these wrappers are always going to be called from outside of the tuning scope
                    if name.startswith("atm") and name.endswith("_work"):
                        continue
                    else:
                        name = name + "_wrapper"
    
                # save
                if name not in t:
                    t[name] = time
    
            for proc_name in t.keys():
                if proc_name not in targeted_subset_runtimes:
                    targeted_subset_runtimes[proc_name] = [t[proc_name]]
                else:
                    targeted_subset_runtimes[proc_name].append(t[proc_name])
        
        except GPTLParseFail as e:
            print(e)

    total = 0
    for proc_name in targeted_subset_runtimes.keys():
        targeted_subset_runtimes[proc_name] = np.mean(targeted_subset_runtimes[proc_name])
        if "::" in proc_name:
            total += targeted_subset_runtimes[proc_name]
        elif proc_name[:proc_name.rfind("_")] in set([n[n.rfind(":") + 1:] for n in t.keys() if "::" in n]):
            total += targeted_subset_runtimes[proc_name]

    with open(os.path.join(working_dir, "gptl_subset_info.pckl"), "wb") as f:
        pickle.dump(targeted_subset_runtimes, f)

    if pprint:
        from pprint import pprint as pretty_print
        pretty_print(sorted(targeted_subset_runtimes.items(), key=lambda item : item[1], reverse=True))

    return total


def gptl_parse_all(working_dir="."):
    procedure_self_runtimes = {}
    module_self_runtimes = {}
    n_timing_files = 0
    for gptl_timing_file in glob(os.path.join(working_dir, "timing.*")):
        
        n_timing_files += 1

        try:
            with open(gptl_timing_file, "r", errors='replace') as f:
                lines = [line.replace(") [0x","__[0x", 1) for line in f.readlines()]

            # find idx for start of subset timing information
            subset_timing_idx = 0
            while subset_timing_idx + 1 < len(lines) and not lines[subset_timing_idx].strip().startswith("Stats for thread 0:"):
                subset_timing_idx += 1
            subset_timing_idx += 1

            # pass over subset timing information to count the number of procedures
            # also, check if the gptl log file mis-printed; if so, we will skip this particular log
            procedure_names = set(["::GPTL_ROOT"])
            i = subset_timing_idx
            line_length_invariant = len(lines[subset_timing_idx])
            while i + 1 < len(lines) and not lines[i + 1].startswith("Overhead sum ="):
                i += 1
                if len(lines[i]) != line_length_invariant:
                    raise GPTLParseFail(f"Parsing of {gptl_timing_file} failed on line {i} while parsing timing info:\n{lines[i]}")

                name = lines[i][2:].split()[0]
                procedure_names.add(get_scoped_name(name))

            # contstruct call matrix
            C_counts = np.zeros((len(procedure_names),len(procedure_names)))
            C_counts[0,0] = 1

            # gather long-name-translation information
            long_name_to_alias_map = {}
            alias_to_long_name_map = {}
            while i + 1 < len(lines) and not lines[i].strip().startswith("thread 0 long name translations"):
                i += 1
            while i + 1 < len(lines) and not lines[i + 1].strip().startswith("Multiple parent info for thread 0:"):
                i += 1
                if " = " in lines[i]: 
                    alias, long_name = (x.strip() for x in lines[i].split("=", maxsplit=1))
                    long_name_to_alias_map[long_name] = alias
                    alias_to_long_name_map[alias] = long_name

            # pass over subset timing information again to map scoped names to call counts 
            procedure_names = ["::GPTL_ROOT"]
            scope_stack = ["::GPTL_ROOT"]
            total_times = np.array([np.nan] * C_counts.shape[0])
            total_times[0] = 0
            i = subset_timing_idx
            line_length_invariant = len(lines[subset_timing_idx])
            while i + 1 < len(lines) and not lines[i + 1].startswith("Overhead sum ="):
                i += 1

                try:
                    callee_name = lines[i][2:].split()[0]
                    call_count = int(float(lines[i][2:].split()[1]))
                    total_time = float(lines[i][2:].split()[3])
                except ValueError:
                    raise GPTLParseFail(f"Parsing of {gptl_timing_file} failed on line {i} while parsing timing info:\n{lines[i]}")

                # check for long-name translation
                if callee_name in long_name_to_alias_map:
                    callee_name = long_name_to_alias_map[callee_name]

                callee_name = get_scoped_name(callee_name)
                if callee_name not in procedure_names:
                    procedure_names.append(callee_name)
                scope_stack.append(callee_name)
                caller_name = scope_stack[-2]

                callee_idx = procedure_names.index(callee_name)
                caller_idx = procedure_names.index(caller_name)

                if np.isnan(total_times[callee_idx]):
                    total_times[callee_idx] = total_time
                elif total_times[callee_idx] != total_time:
                    raise GPTLParseFail(f"Parsing of {gptl_timing_file} failed because of mismatching times for {callee_name}: {total_time} vs {total_times[callee_idx]}")

                multi_parent = lines[i].startswith("*")
                if not multi_parent:
                    C_counts[caller_idx, callee_idx] = call_count

                nesting_level = (len(lines[i+1][2:]) - len(lines[i+1][2:].lstrip())) // 2
                scope_stack = scope_stack[:nesting_level + 1]

            # gather multi-parent information
            while i + 1 < len(lines) and not lines[i + 1].strip().startswith("Multiple parent info for thread 0:"):
                i += 1
            while i + 1 < len(lines) and not lines[i + 1].strip().startswith("Total GPTL memory usage"):
                i += 1
                if lines[i].strip() == "":
                    call_counts = {}
                    while i + 1 < len(lines) and lines[i + 1].strip() != "" and not lines[i + 1].strip().startswith("Total GPTL memory usage"):
                        i += 1

                        try:
                            call_count, proc_name = lines[i].split()
                            call_count = int(float(call_count))
                        except ValueError:
                            raise GPTLParseFail(f"Parsing of {gptl_timing_file} failed on line {i} while parsing multi-parent info:\n{lines[i]}")

                        # check for long-name translation
                        if proc_name in long_name_to_alias_map:
                            proc_name = long_name_to_alias_map[proc_name]

                        proc_name = get_scoped_name(proc_name)

                        # if this isn't the callee entry, save another caller
                        if lines[i + 1].strip() != "":
                            caller_name = proc_name
                            call_counts[caller_name] = call_count

                        # otherwise, save the call counts
                        else:
                            callee_name = proc_name
                            callee_idx = procedure_names.index(callee_name)
                            for caller_name in call_counts.keys():

                                if caller_name == callee_name:
                                    assert(False)
                                else:
                                    try:
                                        caller_idx = procedure_names.index(caller_name)
                                    except:
                                        raise GPTLParseFail(f"Parsing of {gptl_timing_file} failed; procedure {caller_name} in multi-parent info not found in timing info")
                                    C_counts[caller_idx, callee_idx] = call_counts[caller_name]

# - each row of `total_times_stacked` is the total times of all of the procedures
# - `C_counts`
#     - rows are the number of calls to each callee from the caller procedure represented by that row
#         - sums to total number of calls made by each caller procedure
#     - columns are the number of calls made from each caller to the callee procedure represented by that column 
#         - sums to total number of calls to each callee procedure

# - `C_normalized_counts`
#     - columns sum to 1


# - `C_times`
#     - rows give the time spent in each of the callees when called from the procedure represented by that row
#         - PROBLEM: _may sum to more than the total time spent in the caller_
#     - columns give the time spent in the callee procedure represented by that column when called by each caller
#         - _should sum to the total time spent in the callee_

# - `C_times_adjusted`
#     - rows give the _adjusted_ time spent in each of the callees when called from the procedure represented by that row
#         - _sums to less than or equal to the total time spent in the caller_
#     - columns give the time spent in the callee procedure represented by that column when called by each caller
#         - _should sum to the total time spent in the callee_

            C_times_adjusted = np.zeros(C_counts.shape)
            total_times_mod = deepcopy(total_times)

            while True:
                n_zeros = np.count_nonzero(C_times_adjusted)

                C_normalized_counts = np.where(C_counts != 0, (C_counts / np.sum(C_counts, axis=0)), 0)

                total_times_stacked = np.stack([total_times_mod]*C_counts.shape[0])

                C_times = C_normalized_counts * total_times_stacked
                
                C_times_adjusted += np.where(C_times > total_times_stacked.T, total_times_stacked.T, 0)

                if np.count_nonzero(C_times_adjusted) == n_zeros:
                    break
                else:
                    C_counts = np.where(C_times > total_times_stacked.T, 0, C_counts)
                    total_times_mod = total_times - np.sum(C_times_adjusted, axis=0)

            C_times_adjusted += np.where(C_counts != 0, (C_counts / np.sum(C_counts, axis=0)), 0) * total_times_mod
            self_times = total_times - np.sum(C_times_adjusted, axis=1)

            for proc_idx, proc_name in enumerate(procedure_names):
                if proc_name in procedure_self_runtimes:
                    procedure_self_runtimes[proc_name].append(self_times[proc_idx])
                else:
                    procedure_self_runtimes[proc_name] = [self_times[proc_idx]]

        except GPTLParseFail as e:
            print(e)
            n_timing_files -= 1

    overall_runtime = 0
    for proc_name in procedure_self_runtimes.keys():

        procedure_self_runtimes[proc_name] = np.mean(procedure_self_runtimes[proc_name])
        if proc_name != "::GPTL_ROOT":
            overall_runtime += procedure_self_runtimes[proc_name]

        if proc_name.count("::") > 1:
            module_name = proc_name[: 2 + proc_name[2:].find("::")]
        else:
            module_name = proc_name

        if module_name in module_self_runtimes:
            module_self_runtimes[module_name] += procedure_self_runtimes[proc_name]
        else:
            module_self_runtimes[module_name] = procedure_self_runtimes[proc_name]

    print(f"total CPU time: {overall_runtime}")
    print("===================================================================")
    print("module name                                    self           ")
    for module_name in sorted(module_self_runtimes.keys(), key = lambda x: module_self_runtimes[x], reverse=True):
        print(f"{module_name:<47}{module_self_runtimes[module_name]:.2e} ({module_self_runtimes[module_name]/overall_runtime * 100:.2f}%)")


    print("===================================================================")
    print("procedure name                                    self           ")
    for procedure_name in sorted(procedure_self_runtimes.keys(), key = lambda x: procedure_self_runtimes[x], reverse=True):
        print(f"{procedure_name:<47}{procedure_self_runtimes[procedure_name]:.2e} ({procedure_self_runtimes[procedure_name]/overall_runtime * 100:.2f}%)")


def get_scoped_name(name, append_dummy_variable=False):
    '''
    change name format from gptl printout to the scoped format used by PROSE for variables, possibly completed with a dummy variable name
    '''

    scoped_name = name
    if scoped_name.startswith("__"):
        scoped_name = scoped_name[2:]
    scoped_name = "::" + scoped_name
    scoped_name = scoped_name.replace("_MOD_", "::")
    scoped_name = scoped_name.replace("module_mp_", "module_prosemp_") # added as a workaround for mpas' module_mp_ naming convention
    scoped_name = scoped_name.replace("_mp_", "::")
    scoped_name = scoped_name.replace("_prosemp_", "_mp_")    
    if scoped_name.endswith("_"):
        scoped_name = scoped_name[:-1]
    if append_dummy_variable:
        scoped_name =  scoped_name + "::prose"

    return scoped_name


if __name__ == "__main__":
    try:
        print(
            gptl_parse_subset(
                search_patterns=["^::ITPACKV::.+"],
                ignore_patterns=["^::MESSENGER::.+"],
                working_dir=".",
                pprint = True
            )
        )
        # gptl_parse_all()
    except:
        import pdb, traceback, sys
        _, _, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)