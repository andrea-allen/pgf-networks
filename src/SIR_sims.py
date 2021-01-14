import networkx as nx
import numpy as np
import math
from src import event_driven
import time

# TODO have the pgf formalism file only contain "do" not visualize or analyze

# TODO make results easily returnable
# TODO deal with display and pos, which won't work or be solved for large networks

def run():

    # Manipulate-able method for running whatever simulations and plotting we want

    # Assign degree distribution:
    degree_distrb = power_law_degree_distrb()
    degree_distrb = binomial_degree_distb(100)

    # Run two sets of ensembles: one base level with no intervention, one with intervention introduced by specified params
    simulate_intervention_effects(degree_distrb, 'power_law_08_to04_gen3_fast', 2000, 1000,
                                                               0.15, 4, 0.05, .001)

    # Runs 50,000 simulations with a power law degree distribution
    degree_distrb = power_law_degree_distrb()
    # simulate_and_compare_rounds_with_with_without_intervention(degree_distrb, 'power_law_08_to_06_gen3', 50000, 1000, 0.8, 3, 0.6, .001)

    degree_distrb = binomial_degree_distb(1000)

    # Runs 50,000 simulations with binomial degree distribution, saves results
    # simulate_and_compare_rounds_with_with_without_intervention(degree_distrb, 'binomial_02_01_gen3', 50000, 1000, 0.2, 3, 0.1, .001)

    # simulate_and_compare_rounds_with_with_without_intervention(degree_distrb, 'binomial_02_01_gen4', 50000, 1000, 0.2, 4, 0.1, .001)
    print('done')


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
    sim = event_driven.Simulation(1000000, G, beta, gamma, Lambda, Gamma, pos, A)
    sim.run_sim(intervention_gen, beta_interv)
    if current % 5 == 0:
        print('current sim ' + str(current))
        print("--- %s seconds to run simulation---" % (time.time() - start_time))
    results = sim.tabulate_observables(100)
    return results


def power_law_degree_distrb():
    degree_dist = np.zeros(40)
    for k in range(1, len(degree_dist)):
        p_k = (k ** (-2)) * (math.e ** (-k / 5))
        degree_dist[k] = p_k
    return degree_dist


def binomial_degree_distb(N):
    degree_dist = np.zeros(40)
    p = 6 / N
    for k in range(0, len(degree_dist)):
        p_k = (p ** k) * ((1 - p) ** (N - k)) * math.comb(N, k)
        degree_dist[k] = p_k
    return degree_dist


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
