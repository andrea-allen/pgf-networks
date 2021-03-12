import networkx as nx
import numpy as np
import math
from src import event_driven
import time

# Phasing out this file to make way for ensemble.py, which will call the external epintervene package

def simulate_intervention_effects(degree_distrb, base_file_name='sim_results',
                                                               num_sims=10000, num_nodes=1000, init_T=0.8,
                                                               gen_intervene=3, T_intervene=0.4, recover_rate=.001):
    """

    :rtype: void
    """
    # WITH NO intervention
    start_time = time.time()
    size_distrb_per_gen_no_intervention = simulate_ensemble(degree_distrb, num_sims, num_nodes,
                                                                                     -1, 0.0, init_T, recover_rate)
    print('Total time for ensemble was ', time.time() - start_time, ' with N=', num_nodes, ' nodes and ', num_sims,
          ' number of simulations.')
    # TODO have the simulations writing to a file every 100 or so simulations instead of doing it all in memory, as a toggle option
    np.savetxt(base_file_name + '_size_distrb_per_gen_no_interv.txt', size_distrb_per_gen_no_intervention,
               delimiter=',')

    # WITH intervention
    start_time = time.time()
    size_distrb_per_gen_intervention = simulate_ensemble(degree_distrb, num_sims, num_nodes,
                                                                               gen_intervene, T_intervene, init_T,
                                                                               recover_rate)
    print('Total time for ensemble was ', time.time() - start_time, ' with N=', num_nodes, ' nodes and ', num_sims,
          ' number of simulations.')
    np.savetxt(base_file_name + '_size_distrb_per_gen_with_interv.txt', size_distrb_per_gen_intervention, delimiter=',')


def simulate_ensemble(degree_distrb, num_sims=10, N=1000, intervention_gen=-1, intervention_T=0.0, initial_T=0.8,
                      gamma=0.1):

    """

    :rtype: bytearray
    """
    # If intervention_gen is set to -1, no intervention will occur
    # Generates the empirical s-slice for results of proportion of simulations that resulted in exactly s nodes infected at generation g
    # Includes steps for doing so with an intervention step
    beta_init = -(gamma * initial_T) / (initial_T - 1)
    beta_interv = -(gamma * intervention_T) / (intervention_T - 1)
    # s_sizes = np.arange(N)
    outbreak_size_distrb_per_gen_matrix = np.zeros((100, N))
    G, pos = generate_graph(N, degree_distrb)
    A = np.array(nx.adjacency_matrix(G).todense())
    num_nodes_in_net = len(A[0])
    Lambda = np.full((num_nodes_in_net, num_nodes_in_net), beta_init)
    Gamma = np.full(num_nodes_in_net, gamma)
    for i in range(num_sims):
        if i % 500 == 0:
            G, pos = generate_graph(N, degree_distrb)
            A = np.array(nx.adjacency_matrix(G).todense())
            num_nodes_in_net = len(A[0])
            Lambda = np.full((num_nodes_in_net, num_nodes_in_net), beta_init)
            Gamma = np.full(num_nodes_in_net, gamma)
        results = simulate(G, A, pos, beta_init, gamma, Lambda, Gamma, i, intervention_gen, beta_interv)
        for gen in range(len(results)):
            num_total_infctd = int(results[gen])
            try:
                outbreak_size_distrb_per_gen_matrix[gen][num_total_infctd] += 1
            except IndexError:
                print('Index error for g: ', gen, ', gen_s: ', num_total_infctd)
                continue
    # averaging:
    for gen in range(100):
        gen_time_series = outbreak_size_distrb_per_gen_matrix[gen]
        gen_time_series = gen_time_series / num_sims
        outbreak_size_distrb_per_gen_matrix[gen] = gen_time_series
    return outbreak_size_distrb_per_gen_matrix


def simulate(G, A, pos, beta, gamma, Lambda, Gamma, current, intervention_gen=-1, beta_interv=-1.0):
    start_time = time.time()
    # With intervention into the simulation code
    sim = None
    # If intervention is going to be applied then create a UniversalIntervention class
    if intervention_gen > 0:
        # sim = event_driven.UniversalInterventionSim(1000000, G, beta, gamma, Lambda, Gamma, pos, A, intervention_gen, beta_interv)
        # Example of random vaccination:
        sim = event_driven.RandomInterventionSim(1000000, G, beta, gamma, Lambda, Gamma, pos, A, intervention_gen, beta_interv, 0.4)
    else:
        sim = event_driven.Simulation(1000000, G, beta, gamma, Lambda, Gamma, pos, A)
    sim.run_sim()
    if current % 50 == 0:
        print('current sim ' + str(current))
        print("--- %s seconds to run simulation---" % (time.time() - start_time))
    results = sim.tabulate_observables(100)
    # TODO can add a second set of results here that count for the time series of infectives
    return results


def generate_graph(N, deg_dist):
    start_time_1 = time.time()
    number_of_nodes = N * np.array(deg_dist)
    degree_sequence = []
    for i in range(int(math.floor(len(number_of_nodes)))):
        number_with_that_degree = number_of_nodes[i]
        for k in range(int(math.floor(number_with_that_degree))):
            degree_sequence.append(i)
    graphical = nx.is_graphical(degree_sequence)
    # print('Degree sequence is graphical: ', graphical)
    if not graphical:
        degree_sequence.append(1)
    G = nx.configuration_model(degree_sequence)
    # Remove self-loops and parallel edges
    try:
        G.remove_edges_from(nx.selfloop_edges(G))
    except RuntimeError:
        print('No self loops to remove')
    pos = None
    print("--- %s seconds to create graph ---" % (time.time() - start_time_1))
    # start_time_2 = time.time()
    # pos = nx.spring_layout(G)
    # print("--- %s seconds to create pos---" % (time.time() - start_time_2))
    # nx.draw_networkx_nodes(G, pos=pos, with_labels=True)
    # nx.draw_networkx_labels(G, pos=pos, with_labels=True)
    # nx.draw_networkx_edges(G, pos=pos)
    # plt.show()
    # Check:
    # inferred_degree_dist = np.array(nx.degree_histogram(G)) / N
    # print('Inferred equals given degree distribution: ', inferred_degree_dist == deg_dist)
    return G, pos
