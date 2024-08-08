from itertools import product

class BruteForceSearch:

    def __init__(self, search_space):
        self.completed_config_counter = 0
        self.config_iterator = product(*[list(reversed(sorted(kinds))) for kinds in search_space.values()])
        self.var_names = list(search_space.keys())

    def feedback(self, configuration_dict):        
        self.completed_config_counter += 1

    def get_next(self):

        next_config = next(self.config_iterator, None)

        if next_config:
            return { 'config' : {self.var_names[i] : next_config[i] for i in range(len(next_config))} }
        else:
            return None