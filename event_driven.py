import random
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def draw_tau():
    tau = random.expovariate(1.0)
    return tau


def draw_event_class(V_IS, V_I):
    recovery_rate = 0
    infection_rate = 0
    for node in V_I:  # infected nodes
        recovery_rate += node.event_rate
    for edge in V_IS:
        infection_rate += edge.event_rate
    r_start = 0
    r_end = recovery_rate
    i_end = r_end + infection_rate
    random_draw = random.uniform(r_start, i_end)
    if random_draw < r_end:
        return 2  # for recovery event
    elif (random_draw >= r_end) & (random_draw < i_end):
        return 1  # for infection event


def draw_specific_event(max_rate, event_list):
    accepted = False
    random_event = None
    while not accepted:
        random_event = random.choice(event_list)
        accept_rate = random_event.event_rate / max_rate
        random_draw = random.uniform(0, 1)
        if random_draw < accept_rate:
            accepted = True
    return random_event


class Simulation:
    def __init__(self, total_sim_time, G, Lambda, Gamma, pos):
        self.sim_time = total_sim_time
        self.current_sim_time = 0
        self.G = G
        self.A = None
        self.Lambda = Lambda
        self.Gamma = Gamma
        self.V_IS = []
        self.V_I = []
        self.has_been_infected_labels = []
        self.V_R = []
        self.graph_pos = pos
        self.gen_collection = {}

    def intialize(self):
        self.A = np.array(nx.adjacency_matrix(self.G).todense())
        p_zero_idx = random.randint(0, len(self.A[0])-1)
        patient_zero = Node(p_zero_idx, 0, 1, self.Gamma[p_zero_idx])
        self.gen_collection[0] = [p_zero_idx]
        self.V_I.append(patient_zero)
        self.has_been_infected_labels.append(p_zero_idx)
        for j in range(0, len(self.A[0])):
            if self.A[p_zero_idx, j] == 1:
                edge_ij = Edge(patient_zero, Node(j, -1, 0, self.Gamma[j]), self.Lambda[p_zero_idx, j])
                self.V_IS.append(edge_ij)

    def run_sim(self):
        self.intialize()
        while self.current_sim_time < self.sim_time:
            self.visualize_network()
            print('time: ', self.current_sim_time)
            # print('Current IS edges: ')
            # for edge in self.V_IS:
            #     edge.display_info()
            # print('Current Infected nodes: ')
            # for node in self.V_I:
            #     node.display_info()
            # print('Current Recovered nodes: ')
            # for node in self.V_R:
            #     node.display_info()
            tau = draw_tau()
            event_class = draw_event_class(self.V_IS, self.V_I)
            if event_class == 1:
                infection_event = draw_specific_event(np.max(self.Lambda), self.V_IS)
                infection_event.infect()
                self.V_IS.remove(infection_event)
                self.V_I.append(infection_event.j)
                try:
                    self.gen_collection[infection_event.j.gen].append(infection_event.j.i)
                except KeyError:
                    self.gen_collection[infection_event.j.gen] = [infection_event.j.i]
                self.has_been_infected_labels.append(infection_event.j.i)
                self.add_IS_edges(infection_event.j)
            if event_class == 2:
                recovery_event = draw_specific_event(np.max(self.Gamma), self.V_I)
                self.V_I.remove(recovery_event)
                self.remove_IS_edges(recovery_event)
                self.V_R.append(recovery_event)
            self.current_sim_time += tau
            if len(self.V_IS) == 0:
                break

    def add_IS_edges(self, infected_node):
        for j in range(0, len(self.A[infected_node.i])):
            if self.A[infected_node.i][j] == 1:
                if j not in self.has_been_infected_labels:
                    edge_ij = Edge(infected_node, Node(j, -1, 0, self.Gamma[j]), self.Lambda[infected_node.i][j])
                    self.V_IS.append(edge_ij)

    def remove_IS_edges(self, infected_node):
        print('VIS has num edges, ', len(self.V_IS))
        updated_V_IS = []
        for edge in self.V_IS:
            if edge.i.i != infected_node.i:
                updated_V_IS.append(edge)
        self.V_IS = updated_V_IS

    def visualize_network(self):
        val_map = {}
        max_gen = max(self.gen_collection.keys())+1
        for gen in self.gen_collection.keys():
            nodes = self.gen_collection[gen]
            for node in nodes:
                val_map[node] = gen/max_gen

        values = [val_map.get(node, 0) for node in self.G.nodes()]

        nx.draw_networkx_nodes(self.G, pos=self.graph_pos, cmap=plt.get_cmap('inferno'), node_color=values)
        # nx.draw_networkx_nodes(self.G, pos=self.graph_pos, node_color=)
        nx.draw_networkx_labels(self.G, pos=self.graph_pos, with_labels=True)
        nx.draw_networkx_edges(self.G, pos=self.graph_pos)
        plt.show()
        return 0



class Node:
    def __init__(self, i, gen, state, recover_rate):
        self.gen = gen
        self.i = i
        self.state = state
        self.event_rate = recover_rate

    def infect(self):
        self.state = 1

    def recover(self):
        self.state = 2

    def set_gen(self, g):
        self.gen = g

    def display_info(self):
        print('Node index: ', self.i, ' state: ', self.state, ' event_rate: ', self.event_rate, ' gen: ', self.gen)


class Edge:
    def __init__(self, i, j, infect_rate):
        self.i = i #not just an index, this is a whole Node object
        self.j = j
        self.event_rate = infect_rate

    def infect(self):
        self.j.infect()
        self.j.set_gen(self.i.gen + 1)

    def display_info(self):
        print('Edge with event rate: ', self.event_rate, ' nodes:')
        self.i.display_info()
        self.j.display_info()
