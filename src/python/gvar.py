import numpy as np
import networkx as nx

DEFAULT_KIND = 8

class VariableInteractionGraph:

    def __init__(self, G_var_path, sourceTransformers):
        with open(G_var_path, "r") as f:
            lines = f.readlines()

            # load all vertices, already in vertexD order
            self.var_names = []
            for line in lines:
                if "{" in line or "}" in line:
                    continue    # skip first and last line of .dot file
                elif line.find("--") == -1:
                    try:
                        var_name = line.split('"')[1]
                    except IndexError as e:
                        var_name = line.split('=')[1].split("]")[0]

                    assert( var_name == var_name.lower() )
                    self.var_names.append(var_name)

            matrix = np.zeros((len(self.var_names), len(self.var_names)))

            # load all edges
            for line in lines:
                if line.find("--") != -1:
                    source_vertexD = int(line.split('--')[0].strip())
                    target_vertexD = int(line.split('[')[0].split('--')[1].strip())
                    weight = float(line.split("=")[1].split("]")[0].strip().replace('"',''))
                    matrix[source_vertexD][target_vertexD] += weight
                    matrix[target_vertexD][source_vertexD] += weight

        self.sourceTransformers = sourceTransformers

        # construct graph
        try:
            self.G = nx.from_numpy_matrix(matrix)
        except AttributeError:
            self.G = nx.from_numpy_array(matrix)
            
        nx.relabel_nodes(self.G, {i : name for i, name in enumerate(self.var_names)}, copy=False)
        nx.set_node_attributes(self.G, {name : DEFAULT_KIND for name in self.var_names}, name="kind")
        nx.set_node_attributes(self.G, {name : name[:name.rfind("::")] for name in self.var_names}, name="scope")

        # get map of all of the original edge weights for resetting graph
        self.original_edge_weights = nx.get_edge_attributes(self.G, "weight")
        
        # member variables to be set once search starts
        # self.subgraph_nodes = None
        self.original_cost = None


    def get_cost_ratio(self, configuration_obj):

        cost, _ = self.get_cost(configuration_obj["config"])
        return self.original_cost / cost


    def get_cost(self, configuration, custom={}):
        
        # update note attributes with the new kinds
        nx.set_node_attributes(self.G, configuration, name="kind")
            
        target_nodes = list(configuration.keys())
        subgraph_nodes = set(target_nodes)
        for node in target_nodes:
            subgraph_nodes = subgraph_nodes.union(set(self.G.neighbors(node)))

        # calculate and update the edge weights that result
        deltas = self.update_edge_weights(subgraph_nodes, custom)

        # calculate sum of new edge weights
        cost = self.get_edge_weight_sum(subgraph_nodes)

        # reset
        self.reset_graph()

        # for the first configuration, we must set the original cost
        if self.original_cost == None:
            self.original_cost = cost
        
        return cost, deltas


    def update_edge_weights(self, subgraph_nodes, custom={}):

        intra_mixed_change = 0
        intra_low_change = 0
        inter_mixed_change = 0
        inter_low_change = 0

        if custom:
           intra_mixed = custom["intra_mixed"]
           intra_low = custom["intra_low"]
           inter_mixed = custom["inter_mixed"]
           inter_low = custom["inter_low"]
        else:
            intra_mixed = lambda weight : weight * 2
            intra_low = lambda weight : weight / 2
            inter_mixed = lambda weight, dim : abs(2 * weight * (1 + dim))
            inter_low = lambda weight, dim : abs(weight)

        sub_G = self.G.subgraph(subgraph_nodes)
        edge_weights = nx.get_edge_attributes(sub_G, "weight")
        for edge in sub_G.edges():

            scopedVarName0 = edge[0]
            scopedVarName1 = edge[1]
            node1_kind = int(self.G.nodes[scopedVarName0]['kind'])
            node2_kind = int(self.G.nodes[scopedVarName1]['kind'])

            # mixed precision case
            if node1_kind + node2_kind == 12:
                
                # cost equal to original weight * (1 + dim) for interprocedural flow
                # double cost for intraprocedural flow
                if edge_weights[edge] < 0:
                    target_scope0 = self.G.nodes[scopedVarName0]['scope']
                    dim0 = self.sourceTransformers[target_scope0].variable_profile[scopedVarName0]["dim"]

                    target_scope1 = self.G.nodes[scopedVarName1]['scope']
                    dim1 = self.sourceTransformers[target_scope1].variable_profile[scopedVarName1]["dim"]

                    dim = min(dim0,dim1)

                    new_weight = inter_mixed(edge_weights[edge], dim)
                    inter_mixed_change += (abs(new_weight) - abs(edge_weights[edge]))
                    edge_weights[edge] =  new_weight
                else:
                    new_weight = intra_mixed(edge_weights[edge])
                    intra_mixed_change += (new_weight - edge_weights[edge])
                    edge_weights[edge] =  new_weight

            # reduced precision case
            elif node1_kind + node2_kind == 8:
                
                # no change in cost for interprocedural dataflow
                # half cost for intraprocedural dataflow
                if edge_weights[edge] < 0:
                    target_scope0 = self.G.nodes[scopedVarName0]['scope']
                    dim0 = self.sourceTransformers[target_scope0].variable_profile[scopedVarName0]["dim"]

                    target_scope1 = self.G.nodes[scopedVarName1]['scope']
                    dim1 = self.sourceTransformers[target_scope1].variable_profile[scopedVarName1]["dim"]

                    dim = min(dim0,dim1)
                    new_weight = inter_low(edge_weights[edge], dim)
                    inter_low_change += (abs(new_weight) - abs(edge_weights[edge]))
                    edge_weights[edge] =  new_weight

                else:
                    new_weight = intra_low(edge_weights[edge])
                    intra_low_change += (new_weight - edge_weights[edge])
                    edge_weights[edge] =  new_weight

            # high precision case, no change since this is expected to be default
            elif node1_kind + node2_kind == 16:
                pass

            else:
                assert(False)


        nx.set_edge_attributes(self.G, edge_weights, name="weight")

        return {"intra_mixed_change" : intra_mixed_change,
                "intra_low_change" : intra_low_change,
                "inter_mixed_change" : inter_mixed_change,
                "inter_low_change" : inter_low_change,}


    def get_edge_weight_sum(self, subgraph_nodes):        
        sub_G = self.G.subgraph(subgraph_nodes)
        edge_weights = nx.get_edge_attributes(sub_G, "weight")
        return sum([abs(x) for x in edge_weights.values()])


    def reset_graph(self):
        nx.set_node_attributes(self.G, {name : DEFAULT_KIND for name in self.var_names}, name="kind")
        nx.set_edge_attributes(self.G, self.original_edge_weights, name="weight")