import random
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import time
from deprecated import deprecated


# CLASS DEPRECATED: Features moved over to new package epintervene available on TestPyPI
class Node:
    def __init__(self, label, gen, state, recover_rate):
        self.gen = gen
        self.label = label
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

    def equals(self, other_edge):
        # Imperative to use Node class equality here
        if self.i.equals(other_edge.i) and self.j.equals(other_edge.j):
            return True
        else:
            return False


# TODO rename some of the confusing parameters
# TODO rename everything called Lambda to be sum of rates, and actual lambda should be capital Beta, small lambdas? for particular event rates
@deprecated
class Simulation:
    def __init__(self, total_sim_time, G, beta, gamma, Lambda, Gamma, pos, A):
        self.sim_time = total_sim_time
        self.current_sim_time = 0
        self.G = G
        self.beta = beta
        self.gamma = gamma
        self.A = A
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
        self.time_of_intervention = total_sim_time
        self.total_num_timesteps = 0
        self.max_beta = 0
        self.sum_of_rates_vector = []
        self.tau_values = []
        self.time_series = [0]
        self.real_time_srs_infc = []
        self.real_time_srs_rec = []
        self.generational_emergence = {0 : 0}

    def initialize(self):
        start_time = time.time()
        self.max_beta = np.max(self.Lambda)
        # print('starting beta is, ', self.beta)
        N = len(self.A[0])
        p_zero_idx = random.randint(0, N - 1)
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

    def configure(self):
        print('TODO: Add configure options from parent that can be overriden by child')

    # Eventually make classes that extend this base class that can have interventions applied on top of them?
    # def configure_interventions(self, ):
    def run_sim(self):
        self.initialize()
        while self.current_sim_time < self.sim_time:
            # Run one step
            self.single_step()

            self.total_num_timesteps += 1
            if len(self.V_IS) == 0:
                break

    def single_step(self, visualize=False):
        if visualize:
            self.visualize_network()
        sum_of_rates = determine_draw_tau(self.V_IS, self.V_I, self.beta, self.gamma)
        tau = draw_tau(sum_of_rates)
        self.sum_of_rates_vector.append(sum_of_rates)
        self.tau_values.append(tau)
        self.time_series.append(self.time_series[-1] + tau)
        self.real_time_srs_infc.append(len(self.V_I))
        self.real_time_srs_rec.append(len(self.V_R))
        event_class = draw_event_class(self.V_IS, self.V_I)
        if event_class == 1:
            infection_event = draw_specific_event(self.max_beta, self.V_IS)
            infection_event.infect()
            self.V_IS.remove(infection_event)
            self.V_I.append(infection_event.j)

            try:
                self.gen_collection[infection_event.j.gen].append(infection_event.j.label) #maybe move toward storing actual node objects? but also this could get huge. Could also append a vector that tracks what real time the node became infected too
            except KeyError: # Need a better way than KeyError to catch a new generation
                self.gen_collection[infection_event.j.gen] = [infection_event.j.label]
                self.highest_gen += 1
                self.generational_emergence[self.highest_gen] = self.current_sim_time
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

    def add_IS_edges(self, infected_node):
        for j in range(0, len(self.A[infected_node.label])):
            if self.A[infected_node.label][j] == 1:
                candidate_node = Node(j, -1, 0, self.Gamma[j])
                neighbor_node = self.existing_node(candidate_node)
                if neighbor_node.state == 0:
                    edge_ij = Edge(infected_node, neighbor_node, self.Lambda[infected_node.label][j])
                    if not self.edge_list_contains(edge_ij):
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

    def prune_IS_edges(self):
        for edge in self.V_IS:
            try:
                edge_exists_in_network = (self.A[edge.i.label][edge.j.label] == 1)
                if not edge_exists_in_network:
                    self.V_IS.remove(edge)
            except IndexError:
                self.V_IS.remove(edge)  # This happens if the new network no longer contains that node, can remove them

    def edge_list_contains(self, edge):
        for e in self.V_IS:
            if e.equals(edge):
                return True
        return False

    @deprecated
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

    def tabulate_observables(self, max_gens):
        # todo would be good to also have a "real time" vector that gives the option to plot that version of a timeseries
        # do this next in order to compare time intervention vs not, time steps and real time? just real time probably
        gens = max_gens
        total_infected_time_srs = np.zeros(gens)
        s = 1
        m = 1
        s_max = 1
        for gen in range(max_gens):  # {0: 1, 1: 12, 14, 2: 16, 42, ....
            total_infected_time_srs[gen] = s
            try:
                m = len(self.gen_collection[gen + 1])  # num infected in gen g
                s += m
                s_max = s
            except KeyError:
                # make m=0 and s=the last s for the rest of the "time series"
                s = s_max
                m = 0
        return total_infected_time_srs

    def display_info(self):
        print('V_IS edges: ')
        for edge in self.V_IS:
            edge.display_info()
        print('Infected nodes: ')
        for node in self.V_I:
            node.display_info()

    def simtype(self):
        print('I am a generic Simulation class instance')


# This is just random vaccination with 100% of the population getting the vaccination (but also allows for a uniform reduction in transmission probability)
@deprecated
class UniversalInterventionSim(Simulation):
    def __init__(self, total_sim_time, G, beta, gamma, Lambda, Gamma, pos, A, intervention_gen, beta_interv):
        super().__init__(total_sim_time, G, beta, gamma, Lambda, Gamma, pos, A)
        self.intervention_gen = intervention_gen
        self.beta_interv = beta_interv

    def simtype(self):
        print('I am a simulation class of type universal intervention')

    def run_sim(self):
        self.initialize()
        while self.current_sim_time < self.sim_time:
            # intervention if applicable:
            if not self.intervened:
                if self.highest_gen >= self.intervention_gen:
                    self.intervene(self.beta_interv)
                    self.time_of_intervention = self.current_sim_time
                    self.intervened = True
            # Run one step
            self.single_step()

            self.total_num_timesteps += 1
            if len(self.V_IS) == 0:
                break

    # move to universal intervention subclass
    def intervene(self, reduce_current_edges=False):
        # Simplest intervention is just to re-assign Lambda with uniform the new T value
        # in the future, can hand-select which Lambda to change (vaccinating a fraction of the population)
        N = len(self.Lambda[0])
        new_Lambda = np.full((N, N), self.beta_interv)
        self.Lambda = new_Lambda
        self.beta = self.beta_interv
        self.max_beta = np.max(self.Lambda)
        # change event rate for each existing edge pair
        if reduce_current_edges:
            for edge in self.V_IS:
                edge.event_rate = self.Lambda[edge.i.label][edge.j.label]

@deprecated
class AbsoluteTimeNetworkSwitchSim(Simulation):
    def __init__(self, total_sim_time, G, beta, gamma, Lambda, Gamma, pos, A, new_network_matrix, absolute_time_interv):
        super().__init__(total_sim_time, G, beta, gamma, Lambda, Gamma, pos, A)
        self.A_modified = new_network_matrix
        self.intervention_time = absolute_time_interv

    def simtype(self):
        print('I am a simulation class of type Absolute Time Network Intervention')

    def run_sim(self):
        self.initialize()
        while self.current_sim_time < self.sim_time:
            # time-dependent intervention with modified network TODO will need to make this an option for generational modification too
            # todo should it be for absolute time, or time steps?
            if not self.intervened:
                if self.current_sim_time > self.intervention_time:
                    self.intervene()
                    self.time_of_intervention = self.current_sim_time
                    self.intervened = True
            # Run one step
            self.single_step()

            self.total_num_timesteps += 1
            if len(self.V_IS) == 0:
                break

    def intervene(self):
        # todo need to also modify lambda and gamma here
        self.A = self.A_modified
        N = len(self.A[0])
        new_Lambda = np.full((N, N), self.beta)
        self.Lambda = new_Lambda
        new_Gamma = np.full(N, self.gamma)
        self.Gamma = new_Gamma
        self.prune_IS_edges()
        for node in self.V_I:
            self.add_IS_edges(node)
        self.update_IS_edges()
        print('Modifying network')

@deprecated
class RandomInterventionSim(Simulation):
    def __init__(self, total_sim_time, G, beta, gamma, Lambda, Gamma, pos, A, intervention_gen, beta_redux, proportion_reduced):
        super().__init__(total_sim_time, G, beta, gamma, Lambda, Gamma, pos, A)
        self.intervention_gen = intervention_gen
        self.beta_redux = beta_redux
        self.proportion_reduced = proportion_reduced

    def simtype(self):
        print('I am a simulation class of type Random Intervention')

    def run_sim(self):
        self.initialize()
        while self.current_sim_time < self.sim_time:
            if not self.intervened:
                if self.highest_gen >= self.intervention_gen:
                    self.intervene()
                    self.time_of_intervention = self.current_sim_time
                    self.intervened = True
            # Run one step
            self.single_step()

            self.total_num_timesteps += 1
            if len(self.V_IS) == 0:
                break

    def intervene(self, reduce_current_edges=False):
        N = len(self.A[0])
        frac_of_network = self.proportion_reduced * N
        how_many = 1
        if frac_of_network > 1:
            how_many = int(np.round(frac_of_network, 0))
        random_set = np.random.randint(0, N, how_many)
        for node in random_set:
            self.Lambda[node] = np.full(N, self.beta_redux)
            # TODO column as well?
            self.Lambda[:, node] = np.full(N, self.beta_redux).T
        # change event rate for each existing edge pair
        if reduce_current_edges:
            for edge in self.V_IS:
                edge.event_rate = self.Lambda[edge.i.label][edge.j.label]
        print(f'Randomly vaccinating {self.proportion_reduced} percent of nodes from beta value {self.beta} to {self.beta_redux}')


# TODO these should be inside the Simulation class method
@deprecated
def draw_tau(sum_of_rates):
    # print(sum_of_rates)
    tau = np.random.exponential(1 / sum_of_rates)
    return tau

@deprecated
def draw_event_class(V_IS, V_I):  # [list of IS edges] [list of infected nodes]
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


@deprecated
def draw_specific_event(max_rate, event_list):
    accepted = False
    random_event = None
    L = len(event_list)
    while not accepted:
        random_idx = np.random.randint(0, L)
        random_event = event_list[random_idx]
        accept_rate = random_event.event_rate / max_rate  # ex. Edge.event_rate = .2 TODO need to assign max rate when intervention happens
        random_draw = random.uniform(0, 1)
        if random_draw < accept_rate:
            accepted = True
    return random_event

@deprecated
def determine_draw_tau(V_IS, V_I, beta, gamma):
    v_is_count = len(V_IS)
    v_i_count = len(V_I)
    sum_of_rates = v_is_count * beta + v_i_count * gamma
    return sum_of_rates
