import random
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

class Node:
    def __init__(self, i, gen, state, recover_rate):
        self.gen = gen
        self.label = i
        self.state = state
        self.event_rate = recover_rate

    def infect(self):
        self.state = 1

    def recover(self):
        self.state = 2

    def set_gen(self, g):
        self.gen = g

    def display_info(self):
        print('Node index: ', self.label, ' state: ', self.state, ' event_rate: ', self.event_rate, ' gen: ', self.gen)

    def equals(self, node):
        if self.label == node.label:
            return True
        else:
            return False


class Edge:
    def __init__(self, i, j, infect_rate):
        self.i = i  # not just an index, this is a whole Node object
        self.j = j
        self.event_rate = infect_rate

    def infect(self):
        self.j.infect()
        self.j.set_gen(self.i.gen + 1)

    def display_info(self):
        print('Edge with event rate: ', self.event_rate, ' nodes:')
        self.i.display_info()
        self.j.display_info()


class Simulation:
    def __init__(self, total_sim_time, G, Lambda, Gamma, pos):
        self.sim_time = total_sim_time
        self.current_sim_time = 0
        self.G = G
        self.beta = 0
        self.gamma = 0
        self.A = None
        self.Lambda = Lambda
        self.Gamma = Gamma
        self.V_IS = []
        self.V_I = []
        self.has_been_infected_labels = []
        self.V_R = []
        self.graph_pos = pos
        self.gen_collection = {}
        self.nodes = []
        self.highest_gen = 0
        self.intervened = False
        self.total_num_timesteps = 0

    def intialize(self):
        self.A = np.array(nx.adjacency_matrix(self.G).todense())
        self.beta = np.mean(self.Lambda)
        self.gamma = np.mean(self.Gamma)
        print('starting beta is, ', self.beta)
        p_zero_idx = random.randint(0, len(self.A[0]) - 1)
        patient_zero = Node(p_zero_idx, 0, 1, self.Gamma[p_zero_idx])
        self.nodes.append(patient_zero)
        self.gen_collection[0] = [p_zero_idx]
        self.V_I.append(patient_zero)
        self.has_been_infected_labels.append(p_zero_idx)
        for j in range(0, len(self.A[0])):
            if self.A[p_zero_idx, j] == 1:
                neighbor = Node(j, -1, 0, self.Gamma[j])
                self.nodes.append(neighbor)
                edge_ij = Edge(patient_zero, neighbor, self.Lambda[p_zero_idx, j])
                self.V_IS.append(edge_ij)

    def run_sim(self, intervention_gen=-1, beta_interv=0.0, visualize=False):
        self.intialize()
        while self.current_sim_time < self.sim_time:
            if visualize:
                self.visualize_network()
            # print('time: ', self.current_sim_time)
            # print('Current IS edges: ')
            # for edge in self.V_IS:
            #     edge.display_info()
            # print('Current Infected nodes: ')
            # for node in self.V_I:
            #     node.display_info()
            # print('Current Recovered nodes: ')
            # for node in self.V_R:
            #     node.display_info()

            # intervention if applicable:
            if (not self.intervened) & (intervention_gen > 0):
                if self.highest_gen >= intervention_gen:
                    self.intervene(beta_interv)
                    self.intervened = True

            sum_of_rates = determine_draw_tau(self.V_IS, self.V_I, self.beta, self.gamma)
            tau = draw_tau(sum_of_rates)
            # print(self.current_sim_time)
            event_class = draw_event_class(self.V_IS, self.V_I)
            if event_class == 1:
                infection_event = draw_specific_event(np.max(self.Lambda), self.V_IS)
                infection_event.infect()
                self.V_IS.remove(infection_event)
                self.V_I.append(infection_event.j)
                try:
                    self.gen_collection[infection_event.j.gen].append(infection_event.j.label)
                except KeyError:
                    self.gen_collection[infection_event.j.gen] = [infection_event.j.label]
                    self.highest_gen += 1
                self.has_been_infected_labels.append(infection_event.j.label)
                self.update_IS_edges()
                self.add_IS_edges(infection_event.j)
            if event_class == 2:
                recovery_event = draw_specific_event(np.max(self.Gamma), self.V_I)
                self.V_I.remove(recovery_event)
                recovery_event.recover()
                self.update_IS_edges()
                self.V_R.append(recovery_event)
            self.current_sim_time += tau
            self.total_num_timesteps += 1
            if len(self.V_IS) == 0:
                break

    def add_IS_edges(self, infected_node):
        for j in range(0, len(self.A[infected_node.label])):
            if self.A[infected_node.label][j] == 1:
                candidate_node = Node(j, -1, 0, self.Gamma[j])
                neighbor_node = self.existing_node(candidate_node)
                if neighbor_node.state == 0:
                    edge_ij = Edge(infected_node, neighbor_node, self.Lambda[infected_node.label][j])
                    self.V_IS.append(edge_ij)

    def existing_node(self, candidate_node):
        for node in self.nodes:
            if candidate_node.equals(node):
                return node
        self.nodes.append(candidate_node)
        return candidate_node

    def update_IS_edges(self):
        updated_V_IS = []
        for edge in self.V_IS:
            if (edge.i.state == 1) and (edge.j.state == 0):
                updated_V_IS.append(edge)
        self.V_IS = updated_V_IS

    def visualize_network(self):
        val_map = {}
        max_gen = max(self.gen_collection.keys()) + 1
        for gen in self.gen_collection.keys():
            nodes = self.gen_collection[gen]
            for node in nodes:
                val_map[node] = (gen + 3) / max_gen

        values = [val_map.get(node, 0) for node in self.G.nodes()]

        nx.draw_networkx_nodes(self.G, pos=self.graph_pos, cmap=plt.get_cmap('Reds'), node_color=values)
        # nx.draw_networkx_nodes(self.G, pos=self.graph_pos, node_color=)
        nx.draw_networkx_labels(self.G, pos=self.graph_pos, with_labels=True)
        nx.draw_networkx_edges(self.G, pos=self.graph_pos, edge_color='Grey', lw=2)
        plt.show()
        return 0

    def total_infect_over_all_gens(self, max_gens):
        gens = max_gens
        results_vec = np.zeros((2, gens))
        s = 1
        m = 1
        s_max = 1
        for gen in range(max_gens): #{0: 1, 1: 12, 14, 2: 16, 42, ....
            results_vec[0][gen] = m
            results_vec[1][gen] = s
            try:
                m = len(self.gen_collection[gen + 1])  # num infected in gen g
                s += m
                s_max = s
            except KeyError:
                # make m=0 and s=the last s for the rest of the "time series"
                s = s_max
                m = 0
        return results_vec

    def intervene(self, beta_interv):
        print('intervention')
        # Simplest intervention is just to re-assign Lambda with uniform the new T value
        # in the future, can hand-select which Lambda to change (vaccinating a fraction of the population)
        N = len(self.Lambda[0])
        new_Lambda = np.full((N, N), beta_interv)
        # for i in range(N):
        #     for j in range(N):
        #         new_Lambda[i][j] = beta_interv
        self.Lambda = new_Lambda
        self.beta = beta_interv
        print('new beta', self.beta)
        # change event rate for each existing edge pair
        for edge in self.V_IS:
            edge.event_rate = self.Lambda[edge.i.label][edge.j.label]

def draw_tau(sum_of_rates):
    tau = np.random.exponential(1/sum_of_rates)
    return tau


def draw_event_class(V_IS, V_I): #[list of IS edges] [list of infected nodes]
    recovery_rate = 0
    infection_rate = 0
    for node in V_I:
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
        accept_rate = random_event.event_rate / max_rate #ex. Edge.event_rate = .2
        random_draw = random.uniform(0, 1)
        if random_draw < accept_rate:
            accepted = True
    return random_event

def determine_draw_tau(V_IS, V_I, beta, gamma):
    v_is_count = len(V_IS)
    v_i_count = len(V_I)
    sum_of_rates = v_is_count*beta + v_i_count*gamma
    return sum_of_rates

