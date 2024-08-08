import os
import subprocess
from configparser import ConfigParser, ExtendedInterpolation


class SetupParser(object):
    '''
    Setup File Parser for Prose
    '''

    def __init__(self, filepath, working_dir):

        # initialize data members:
        self._data = None
        self._env = {}
        self.working_dir = working_dir

        # first, parse the setup.ini file and generate self._data of type ConfigParser
        self._parse(filepath)

        # now, run all the env commands and keep a record of env's for each section.
        for section in self._data:
            self._get_env(section)

    @property
    def env(self):
        return self._env

    def __getitem__(self, key):
        if key == 'env':
            return self._env[key]
        else:
            return self._data[key]

    def _parse(self, filepath):
        ''' Given a setup.ini file path, parses the file and stores the data in self._data'''

        self._data = ConfigParser(interpolation=ExtendedInterpolation())
        assert os.path.exists(filepath), "Cannot find setup.ini file at "+filepath
        self._data.read(filepath)

        # process machine-specific section
        assert 'machine' in self._data, "'machine' section is missing in setup file"
        assert 'project_root' in self._data['machine'], "'project_root' entry is missing in 'machine' section of setup file"
        
        self._data['machine']['project_root'] = os.path.abspath(self._data['machine']['project_root'])
        if 'src_search_paths' in self._data['machine']:
            self._data['machine']['src_search_paths'] = "|".join([os.path.abspath(path) for path in self._data['machine']['src_search_paths'].split()])
        else:
            self._data['machine']['src_search_paths'] = self._data['machine']['project_root']

        # process target section
        assert 'target' in self._data, "'target' section is missing in setup file"
        assert 'src_files' in self._data['target'], "'src_files' entry is missing in 'target' section of setup file"
        assert 'search_patterns' in self._data['target'], "'search_patterns' entry is missing in 'target' section of setup file"

        self._data['target']['src_files'] = "|".join([os.path.relpath(path, start=self.working_dir) for path in self._data['target']['src_files'].split()])

        self._data['target']['search_patterns'] = self._data['target']['search_patterns'].strip().lower().replace("\n", "|")
        if 'ignore_patterns' in self._data['target']:
            self._data['target']['ignore_patterns'] = self._data['target']['ignore_patterns'].strip().lower().replace("\n", "|")
        else:
            self._data['target']['ignore_patterns'] = ""

        if 'additional_plugin_flags' in self._data['target']:
            # because the source will be lowered prior to parsing by the plugin, lower the macro define compile flags
            additional_plugin_flags = []
            for plugin_flag in [x.strip() for x in self._data['target']['additional_plugin_flags'].split()]:
                if plugin_flag.startswith("-D"):
                    plugin_flag = plugin_flag[:2] + plugin_flag[2:].lower()
                additional_plugin_flags.append(plugin_flag)
            self._data['target']['additional_plugin_flags'] = " ".join(additional_plugin_flags)
        else:
            self._data['target']['additional_plugin_flags'] = ""

        # process build section
        assert 'build' in self._data, "'build' section is missing in setup file"
        assert 'cmd' in self._data['build'], "'cmd' entry is missing in 'build' section of setup file"

        if 'working_dir' in self._data['build']:
            self._data['build']['working_dir'] = os.path.relpath(self._data['build']['working_dir'], start=self.working_dir)
        else:
            self._data['build']['working_dir'] = os.path.relpath(self._data['machine']['project_root'], start=self.working_dir)
        if 'partial_build_cmd' not in self._data['build']:
            self._data['build']['partial_build_cmd'] = self._data['build']['cmd']

        # process run section
        assert 'run' in self._data, "'run' section is missing in setup file"
        assert 'cmd' in self._data['run'], "'cmd' entry is missing in 'run' section of setup file"

        if 'execution_filtering' in self._data['run']:
            assert self._data['run']['execution_filtering'].lower() in ['true', 'false'], "'execution_filtering' entry is a boolean; specify 'True' or 'False')"
            self._data['run']['execution_filtering'] = str(self._data['run']['execution_filtering'].lower() == 'true').lower()
        else:
            self._data['run']['execution_filtering'] = 'false'
        if 'working_dir' in self._data['run']:
            self._data['run']['working_dir'] = os.path.relpath(self._data['run']['working_dir'], start=self.working_dir)
        else:
            self._data['run']['working_dir'] = os.path.relpath(self._data['machine']['project_root'], start=self.working_dir)
        if not 'timeout' in self._data['run']:
            self._data['run']['timeout'] = '0'

        # process eval section
        assert 'eval' in self._data, "'eval' section is missing in setup file"
        assert 'cmd' in self._data['eval'], "'cmd' entry is missing in 'eval' section of setup file"
        assert 'pass_log_path' in self._data['eval'], "'pass_log_path' entry is missing in 'eval' section of setup file"

        if 'working_dir' in self._data['eval']:
            self._data['eval']['working_dir'] =  os.path.relpath(self._data['eval']['working_dir'], start=self.working_dir)
        else:
            self._data['eval']['working_dir'] = os.path.relpath(self._data['machine']['project_root'], start=self.working_dir)
        if 'cost_threshold' not in self._data['eval']:
            self._data['eval']['cost_threshold'] = "-1.0"

        # process Derecho section if present
        if 'Derecho' in self._data:
            assert 'env_script' in self._data['Derecho'], "'set_env' section is missing in the 'Derecho' section of setup file"
            if 'copy_ignore' in self._data['Derecho']:
                self._data['Derecho']['copy_ignore'] = "|".join([name.strip() for name in self._data['Derecho']['copy_ignore'].strip().split("\n")])
            else:
                self._data['Derecho']['copy_ignore'] = ""


    def _get_env(self, section):
        ''' Given a setup.ini section, runs the env_cmd and stores resulting environment in self.env'''

        if not 'env_cmd' in self._data[section]:
            self._env[section] = None
        else:
            assert section in ['target', 'build', 'run', 'eval'], 'Cannot get env for section '+section

            if self._data[section]['env_cmd'].endswith('.sh') and \
                not (self._data[section]['env_cmd'].startswith('source') or \
                     self._data[section]['env_cmd'].startswith('. ') ):
                raise RuntimeError(section+" env_cmd is a shell script, but is not preceded with "+
                                  "'source' or '.', so will likely be ineffective")

            env = {}
            run = subprocess.run(self._data[section]['env_cmd']+'; env', shell=True, check=True,
                                capture_output=True, encoding='utf-8')
            for line in run.stdout.splitlines():
                if '=' in line[1:-1]:
                    key = line.split('=')[0]
                    val = line.split('=')[1]
                    if key.startswith('BASH_FUNC_'):
                        continue
                    env[key]=val
            self._env[section] = env
