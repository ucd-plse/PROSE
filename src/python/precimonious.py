import subprocess
import pickle
import os
from copy import deepcopy
from math import ceil

class PrecimoniousSearch:

    @staticmethod
    def load(path="prose_workspace/__PrecimoniousSearch.pckl"):
        with open(path, "rb") as f:
            return pickle.load(f)

    def __init__(self, search_space):
        self.available_kinds = {}
        self.config_template = {}
        self.config_queue = []
        self.delta_divisions = 1
        self.cost_threshold = -1
        self.completed_config_counter = 0
        self.current_best_configuration = {}
        self.improvement_flag = False

        for var_name, possible_kinds in search_space.items():

            # initialize each variable with its highest kind
            self.config_template[var_name] = possible_kinds[-1]

            # remember all possible kinds associated with each variable
            self.available_kinds[var_name] = possible_kinds

        # queue up the first two configurations (all vars at their highest kind, all vars at their second highest)
        self.config_queue.append({
            "config" : deepcopy(self.config_template),
            "delta" : list(self.available_kinds.keys()),
            "inv" : False
        })

        self._to_second_highest_precision(self.available_kinds.keys())
        self.config_queue.append({
            "config" : deepcopy(self.config_template),
            "delta" : [],
            "inv" : True
        })


    def save(self):
        with open("prose_workspace/__PrecimoniousSearch.pckl", "wb") as f:
            pickle.dump(self, f)


    def _to_highest_precision(self, target_variables):
        for var_name in target_variables:
            if len(self.available_kinds[var_name]) > 0:
                self.config_template[var_name] = self.available_kinds[var_name][-1]


    def _to_second_highest_precision(self, target_variables):
        for var_name in target_variables:
            if len(self.available_kinds[var_name]) > 1:
                self.config_template[var_name] = self.available_kinds[var_name][-2]


    # if any of the variables still have different precision levels
    # that have yet to be explored, we are not done
    def _done(self):
        return len(self._get_remaining_variable_names()) == 0


    def set_cost_threshold(self, time):
        self.cost_threshold = float(time)


    def feedback(self, configuration_dict):
        
        configuration_dict["configuration_number"] = self.completed_config_counter

        if self.completed_config_counter == 0:
            self.current_best_configuration = configuration_dict
            if self.cost_threshold == -1:
                self.cost_threshold = configuration_dict["cost"]
            
        # these conditionals select for the configuration with the greatest number of low-precision variables that
        # satisfies the runtime threshold. In the case of a tie, the lower-cost configuration is selected
        elif configuration_dict["cost"] <= self.cost_threshold:
            if len(configuration_dict["delta"]) < len(self.current_best_configuration["delta"]):
                self.current_best_configuration = configuration_dict
                self.improvement_flag = True
            elif len(configuration_dict["delta"]) == len(self.current_best_configuration["delta"]):
                if configuration_dict["cost"] < self.current_best_configuration["cost"]:
                    self.current_best_configuration = configuration_dict
                    self.improvement_flag = True

        self.completed_config_counter += 1
        self.save()


    def _generate_next_batch(self):

        log_path = "prose_logs/{:0>4}/".format(self.completed_config_counter - 1)

        subprocess.run(["touch", os.path.join(log_path, "FLAG_----------------------------BATCH_END")])

        if self.improvement_flag:
            subprocess.run(["touch", os.path.join(log_path, "FLAG_IMPROVEMENT_DISCOVERED_IN_CONFIG_{:0>4}".format(self.current_best_configuration["configuration_number"]))])
            self.improvement_flag = False
            
            # if the config is a delta inv, reduce divisions by one;
            # otherwise, reset to divisions to two
            if self.current_best_configuration["inv"]:
                self.delta_divisions -= 1
            else:
                self.delta_divisions = 2

        # otherwise, no new improvement was discovered in the last batch so we increase the granularity
        else:
            subprocess.run(["touch", os.path.join(log_path, "FLAG_NO_IMPROVEMENT")])
            self.delta_divisions *= 2

        # if we cannot divide the current delta any further, it is minimal.
        if self.delta_divisions == 0 or self.delta_divisions > len(self.current_best_configuration["delta"]):

            for var_name in self._get_remaining_variable_names():
                
                # remove any variables whose kind was unable to be lowered
                if self.current_best_configuration["config"][var_name] == self.available_kinds[var_name][-1]:
                    del(self.available_kinds[var_name][:])

                # all remaining variables could have their kind lowered to the second
                # highest. Therefore, we remove the highest kind from the kinds available to them
                if len(self.available_kinds[var_name]) > 0:
                    self.available_kinds[var_name].pop(-1)

            # reset the delta set to all variables still to be
            # searched (note that if we reach this point, the search
            # will only continue if there are variables in the search
            # space with more than two precision levels)
            self.current_best_configuration["delta"] = self._get_remaining_variable_names()
            self.delta_divisions = 2

            subprocess.run(["touch", os.path.join(log_path, "FLAG_RESET_DELTA")])

        subprocess.run(["touch", os.path.join(log_path, "FLAG_DELTA_{}_DIV_{}".format(len(self.current_best_configuration['delta']), self.delta_divisions))])

        if self._done():
            self.config_queue = None
        else:
            self._create_deltas()


    def get_next_batch(self):

        # run config 000 separately
        if self.completed_config_counter == 0:
            return [self.config_queue.pop(0)]
        elif self.config_queue == []:
            self._generate_next_batch()

        ret = self.config_queue
        self.config_queue = []
        return ret


    def get_next(self):

        # if there are no configurations left, create a new batch
        if self.config_queue == []:
            self._generate_next_batch()
                
        # pop the first one
        if self.config_queue:
            return self.config_queue.pop(0)
        else:
            return None


    # return variables that still have multiple prevision levels that
    # have yet to be tested
    def _get_remaining_variable_names(self):
        return [var_name for var_name in list(self.available_kinds.keys()) if len(self.available_kinds[var_name]) > 1]


    def _create_deltas(self):
        old_delta = self.current_best_configuration["delta"]

        div_size = int(ceil(len(old_delta)/self.delta_divisions))

        unique_deltas = set()

        for div_start in range(0, len(old_delta), div_size):
            delta_set = []
            delta_inv_set = []
            for j in range(0, len(old_delta)):
                if j >= div_start and j < div_start + div_size:
                    delta_set.append(old_delta[j])
                else:
                    delta_inv_set.append(old_delta[j])

            if tuple(delta_set) not in unique_deltas:
                unique_deltas.add(tuple(delta_set))
                self._to_second_highest_precision(self._get_remaining_variable_names())
                self._to_highest_precision(delta_set)
                self.config_queue.append({
                    "config" : deepcopy(self.config_template),
                    "delta" : delta_set,
                    "inv" : False
                })

            if tuple(delta_inv_set) not in unique_deltas:
                unique_deltas.add(tuple(delta_inv_set))
                self._to_second_highest_precision(self._get_remaining_variable_names())
                self._to_highest_precision(delta_inv_set)
                self.config_queue.append({
                    "config" : deepcopy(self.config_template),
                    "delta" : delta_inv_set,
                    "inv" : True
                })
