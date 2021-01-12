import networkx as nx
import numpy as np
import math
from src import event_driven
import matplotlib.pyplot as plt
import time

def run():
    # Manipulate-able method for running whatever simulations and plotting we want
    # read_back_data()
    # plt.show()

    # Sims with a power law degree distribution:
    degree_distrb = power_law_degree_distrb()
    #done TODO still need to find out why always slightly less than w the intervention
    #done TODO add toggle for the above logic that should be a feature of the model
    #done TODO need to bring down initialization time from 1.5 secs to much quicker, look into sparse matrix lookups
    #TODO make results easily returnable
    #TODO deal with display and pos, which won't work or be solved for large networks
    #done TODO index issue higher for higher generations?
    simulate_and_compare_rounds_with_with_without_intervention(degree_distrb, 'power_law_08_to04_gen3_fast', 200, 10000, 0.8, 4, 0.0, .001)

    # 1.5 secs to initialize a 10000 node simulation
    # 17 secds for 10000 nodes that "take off"
    # make a distribution of simulation times?

    # Runs 50,000 simulations with a power law degree distribution
    degree_distrb = power_law_degree_distrb()
    # simulate_and_compare_rounds_with_with_without_intervention(degree_distrb, 'power_law_08_to_06_gen3', 50000, 1000, 0.8, 3, 0.6, .001)


    degree_distrb = binomial_degree_distb(1000)

    # Runs 50,000 simulations with binomial degree distribution, saves results
    # simulate_and_compare_rounds_with_with_without_intervention(degree_distrb, 'binomial_02_01_gen3', 50000, 1000, 0.2, 3, 0.1, .001)

    # simulate_and_compare_rounds_with_with_without_intervention(degree_distrb, 'binomial_02_01_gen4', 50000, 1000, 0.2, 4, 0.1, .001)
    print('done')


def simulate_and_compare_rounds_with_with_without_intervention(degree_distrb, base_file_name='sim_results', num_sims=10000, num_nodes=1000, init_T=0.8, gen_intervene=3, T_intervene=0.4, recover_rate=.001):
    # Comment if no need to run or save results:
    # WITH NO intervention
    start_time = time.time()
    s_sizes_no_intervention, size_distrb_per_gen_no_intervention = outbreak_size_distrb_per_gen_with_intervention(degree_distrb, num_sims, num_nodes, -1, 0.0, init_T, recover_rate)
    print('Total time for ensemble was ', time.time() - start_time, ' with N=', num_nodes, ' nodes and ', num_sims, ' number of simulations.')
    np.savetxt(base_file_name+'_size_distrb_per_gen_no_interv.txt', size_distrb_per_gen_no_intervention, delimiter=',')
    # WITH intervention
    start_time = time.time()
    #REAL ONE
    #when initi_T=intervene_T, you get "matching" curves. something's going wrong with the intervention itself possibly, why we need tests
    # s_sizes_intervention, size_distrb_per_gen_intervention = outbreak_size_distrb_per_gen_with_intervention(degree_distrb, num_sims, num_nodes, gen_intervene, T_intervene, init_T, recover_rate)
    s_sizes_intervention, size_distrb_per_gen_intervention = outbreak_size_distrb_per_gen_with_intervention(degree_distrb, num_sims, num_nodes, gen_intervene, T_intervene, init_T, recover_rate)
    # TEST with also no intervention
    # s_sizes_intervention, size_distrb_per_gen_intervention = outbreak_size_distrb_per_gen(degree_distrb, num_sims, num_nodes, init_T, recover_rate)
    print('Total time for ensemble was ', time.time() - start_time, ' with N=', num_nodes, ' nodes and ', num_sims, ' number of simulations.')
    np.savetxt(base_file_name+'_size_distrb_per_gen_with_interv.txt', size_distrb_per_gen_intervention, delimiter=',')

    # Plotting results against one another
    color_key = {1: 'blue', 2: 'red', 3: 'orange', 5: 'black', 7: 'pink', 4:'teal', 20:'magenta'}
    # for gen in [2, 6, 11, 18, 40]:
    for gen in [2, 3, 4, 7, 20]:
        plt.plot(s_sizes_no_intervention[:100], size_distrb_per_gen_no_intervention[gen][:100], label='$g=$' + str(gen), color=color_key[gen], alpha=0.75, ls='-', lw=1)
    plt.legend(loc='upper right')
    plt.xlabel('$s$')
    plt.ylabel('$p_s^g$')
    plt.semilogy()
    # plt.savefig('p_s_g_distribution_no_intervention.png')
    # plt.show()

    # for gen in [2, 6, 11, 18, 40]:
    for gen in [2, 3, 4, 7, 20]:
        plt.plot(s_sizes_intervention[:100], size_distrb_per_gen_intervention[gen][:100], label='$int g=$' + str(gen), color=color_key[gen], alpha=0.75, ls='--', lw=1)
    plt.legend(loc='upper right')
    plt.xlabel('$s$')
    plt.ylabel('$p_s^g$')
    plt.semilogy()
    plt.savefig('p_s_g_distribution_intervention.png')
    plt.show()

    #TODO add an inset plot with total number of infections?

def read_back_data():
    # Manipulatable method for reading back data and plotting desired results

    # data = np.loadtxt('../pgf-nets-data/size_distrb_per_gen_no_int_g3_full.txt', delimiter=',')
    data = np.loadtxt('binomial_02_01_gen3_size_distrb_per_gen_no_interv.txt', delimiter=',')
    # data_int = np.loadtxt('../pgf-nets-data/size_distrb_per_gen_int_g3_full.txt', delimiter=',')
    data_int = np.loadtxt('binomial_02_01_gen3_size_distrb_per_gen_with_interv.txt', delimiter=',')
    color_key = {2: 'blue', 6: 'red', 11: 'orange', 18:'black'}
    for gen in [2, 6, 11]:
        time_series = data[gen][2:200]
        time_series_int = data_int[gen][2:200]
        plt.plot(time_series, label='$g=$' + str(gen), color=color_key[gen], alpha=0.5, lw=1)
        plt.plot(time_series_int, color=color_key[gen], alpha=0.75, ls='--', lw=1)
    plt.legend(loc='upper right')
    plt.xlabel('number infected $s$ at generation $g$')
    plt.ylabel('$p_s^g$')
    plt.title('Effects of simulations with intervention from $T=.2$ to $T=.1$ at $g=3$ on Binomial Degree Distribution Network')
    plt.semilogy()
    plt.ylim(.0001, .1)
    # plt.title('Created from saved data')
    # plt.savefig('p_s_g_distribution_intervention.png')
    plt.show()
    return data


def power_law_degree_distrb():
    degree_dist = np.zeros(40)
    for k in range(1, len(degree_dist)):
        p_k = (k ** (-2)) * (math.e ** (-k / 5))
        degree_dist[k] = p_k
    return degree_dist

def binomial_degree_distb(N):
    degree_dist = np.zeros(40)
    p = 6/N
    for k in range(0, len(degree_dist)):
        p_k = (p**k)*((1-p)**(N-k))*math.comb(N, k)
        degree_dist[k] = p_k
    return degree_dist

def outbreak_size_distrb_per_gen_with_intervention(degree_distrb, num_sims=10, N=1000, intervention_gen=-1, intervention_T=0.0, initial_T=0.8, gamma=0.1):
    # Generates the empirical s-slice for results of proportion of simulations that resulted in exactly s nodes infected at generation g
    # Includes steps for doing so with an intervention step
    beta_init = -(gamma * initial_T) / (initial_T - 1)
    beta_interv = -(gamma*intervention_T)/(intervention_T-1)
    s_sizes = np.arange(N)
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
        for g in range(len(results[0])):
            gen_s = int(results[1][g])
            try:
                outbreak_size_distrb_per_gen_matrix[g][gen_s] += 1
            except IndexError:
                print('Index error for g: ', g, ', gen_s: ', gen_s)
                continue
    # averaging:
    for gen in range(100):
        gen_time_series = outbreak_size_distrb_per_gen_matrix[gen]
        gen_time_series = gen_time_series / num_sims
        outbreak_size_distrb_per_gen_matrix[gen] = gen_time_series
    return s_sizes, outbreak_size_distrb_per_gen_matrix


def simulate(G, A, pos, beta, gamma, Lambda, Gamma, current, intervention_gen=-1, beta_interv=-1.0):
    start_time = time.time()
    # With intervention into the simulation code
    sim = event_driven.Simulation(1000000, G, beta, gamma, Lambda, Gamma, pos, A)
    sim.run_sim(intervention_gen, beta_interv)
    if current % 500 == 0:
        print('current sim ' + str(current))
        print("--- %s seconds to run simulation---" % (time.time() - start_time))
    results = sim.total_infect_over_all_gens(100)
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
