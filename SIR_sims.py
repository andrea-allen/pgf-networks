import networkx as nx
import numpy as np
import math
import event_driven
import matplotlib.pyplot as plt

def run():
    # Manipulate-able method for running whatever simulations and plotting we want
    read_back_data()
    # plt.show()
    # Original degree distribution:
    degree_distrb = power_law_degree_distrb()
    simulate_and_compare_rounds_with_with_without_intervention(degree_distrb)
    print('done')
    degree_distrb = power_law_degree_distrb()
    s_sizes, size_distrb_per_gen = outbreak_size_distrb_per_gen(degree_distrb, 10000, 1000)
    np.savetxt('size_distrb_per_gen.txt', size_distrb_per_gen, delimiter=',')
    for gen in [2, 6, 11, 30]:
        plt.plot(s_sizes[1:350], size_distrb_per_gen[gen][1:350], label='$g=$' + str(gen))
    plt.legend(loc='upper right')
    plt.xlabel('$s$')
    plt.ylabel('$p_s^g$')
    plt.semilogy()
    plt.savefig('p_s_g_distribution.png')
    plt.show()


def simulate_and_compare_rounds_with_with_without_intervention(degree_distrb, num_sims=10000, num_nodes=1000, init_T=0.8, gen_intervene=3, T_intervene=0.4, recover_rate=.001):
    # Comment if no need to run or save results:
    # WITH NO intervention
    s_sizes_no_intervention, size_distrb_per_gen_no_intervention = outbreak_size_distrb_per_gen(degree_distrb, num_sims, num_nodes, init_T, recover_rate)
    np.savetxt('size_distrb_per_gen_no_int_full.txt', size_distrb_per_gen_no_intervention, delimiter=',')
    # WITH intervention
    s_sizes_intervention, size_distrb_per_gen_intervention = outbreak_size_distrb_per_gen_with_intervention(degree_distrb, num_sims, num_nodes, gen_intervene, T_intervene, init_T, recover_rate)
    np.savetxt('size_distrb_per_gen_int_full.txt', size_distrb_per_gen_intervention, delimiter=',')

    # Plotting results against one another
    for gen in [2, 6, 11]:
        plt.plot(s_sizes_no_intervention[2:350], size_distrb_per_gen_no_intervention[gen][2:350], label='$g=$' + str(gen))
    plt.legend(loc='upper right')
    plt.xlabel('$s$')
    plt.ylabel('$p_s^g$')
    plt.semilogy()
    # plt.savefig('p_s_g_distribution_no_intervention.png')
    # plt.show()

    for gen in [2, 6, 11, 18]:
        plt.plot(s_sizes_intervention[2:350], size_distrb_per_gen_intervention[gen][2:350], label='$int g=$' + str(gen))
    plt.legend(loc='upper right')
    plt.xlabel('$s$')
    plt.ylabel('$p_s^g$')
    plt.semilogy()
    plt.savefig('p_s_g_distribution_intervention.png')
    plt.show()

def read_back_data():
    # Manipulatable method for reading back data and plotting desired results

    data = np.loadtxt('../pgf-nets-data/size_distrb_per_gen_no_int_g3_full.txt', delimiter=',')
    data_int = np.loadtxt('../pgf-nets-data/size_distrb_per_gen_int_g3_full.txt', delimiter=',')
    color_key = {2: 'blue', 6: 'red', 11: 'orange', 18:'black'}
    for gen in [2, 6, 11]:
        time_series = data[gen][2:200]
        time_series_int = data_int[gen][2:200]
        plt.plot(time_series, label='$g=$' + str(gen), color=color_key[gen], alpha=0.5, lw=1)
        plt.plot(time_series_int, color=color_key[gen], alpha=0.75, ls='--', lw=1)
    plt.legend(loc='upper right')
    plt.xlabel('number infected $s$ at generation $g$')
    plt.ylabel('$p_s^g$')
    plt.title('Effects of simulations with intervention from $T=.8$ to $T=.4$ at $g=3$')
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


def outbreak_size_distrb_per_gen(degree_distrb, num_sims=10, N=1000, T=0.8, gamma=0.1):
    # Generates the empirical s-slice for results of proportion of simulations that resulted in exactly s nodes infected at generation g
    beta = -(gamma*T)/(T-1)
    s_sizes = np.arange(N)
    outbreak_size_distrb_per_gen_matrix = np.zeros((100, N))
    G, pos = generate_graph(N, degree_distrb)
    N = len(G.nodes())
    Lambda = np.zeros((N, N))
    Gamma = np.zeros(N)
    # Construct recovery values and transmissibility matrix
    for n in range(N):
        Gamma[n] = gamma
        for j in range(N):
            Lambda[n][j] = beta
    # Run simulations, new graph every 500 simulations
    for i in range(num_sims):
        if i % 500 == 0:
            G, pos = generate_graph(N, degree_distrb)
            N_resized = len(G.nodes())
            Lambda = np.zeros((N_resized, N_resized))
            Gamma = np.zeros(N_resized)
            for n in range(N_resized):
                Gamma[n] = gamma
                for j in range(N_resized):
                    Lambda[n][j] = beta
        sm_matrix = simulate_noel(G, pos, Lambda, Gamma, i) #results for time series of s and m nodes infected
        for g in range(len(sm_matrix[0])):
            gen_s = int(sm_matrix[1][g])
            try:
                outbreak_size_distrb_per_gen_matrix[g][gen_s] += 1 #Add to mass of that result
            except IndexError:
                print('Index error for g: ', g, ', gen_s: ', gen_s)
                continue
    # averaging:
    for gen in range(100):
        gen_time_series = outbreak_size_distrb_per_gen_matrix[gen]
        gen_time_series = gen_time_series / num_sims
        outbreak_size_distrb_per_gen_matrix[gen] = gen_time_series
    return s_sizes, outbreak_size_distrb_per_gen_matrix

def outbreak_size_distrb_per_gen_with_intervention(degree_distrb, num_sims=10, N=1000, intervention_gen=-1, intervention_T=0.0, initial_T=0.8, gamma=0.1):
    # Generates the empirical s-slice for results of proportion of simulations that resulted in exactly s nodes infected at generation g
    # Includes steps for doing so with an intervention step
    beta_init = -(gamma * initial_T) / (initial_T - 1)
    beta_interv = -(gamma*intervention_T)/(intervention_T-1)
    s_sizes = np.arange(N)
    outbreak_size_distrb_per_gen_matrix = np.zeros((100, N))
    G, pos = generate_graph(N, degree_distrb)
    N = len(G.nodes())
    Lambda = np.zeros((N, N))
    Gamma = np.zeros(N)
    for n in range(N):
        Gamma[n] = gamma
        for j in range(N):
            Lambda[n][j] = beta_init
    for i in range(num_sims):
        if i % 500 == 0:
            G, pos = generate_graph(N, degree_distrb)
            N_resized = len(G.nodes())
            Lambda = np.zeros((N_resized, N_resized))
            Gamma = np.zeros(N_resized)
            for n in range(N_resized):
                Gamma[n] = gamma
                for j in range(N_resized):
                    Lambda[n][j] = beta_init
        sm_matrix = simulate_noel(G, pos, Lambda, Gamma, i, intervention_gen, beta_interv)
        for g in range(len(sm_matrix[0])):
            gen_s = int(sm_matrix[1][g])
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


def simulate():
    # A sample function for event driven simulation use
    degree_dist = [0, .20, .20, .15, .05, .35, .02, .02, .01]
    N = 30
    G, pos = generate_graph(N, degree_dist)
    Lambda = np.zeros((N, N))
    Gamma = np.zeros(N)
    for i in range(N):
        Gamma[i] = .05
        for j in range(N):
            Lambda[i][j] = .5
    sim = event_driven.Simulation(N * 5, G, Lambda, Gamma, pos)
    sim.run_sim()
    sm_matrix = sim.generate_matrix_gen(20)
    return sm_matrix


def simulate_noel(G, pos, Lambda, Gamma, current, intervention_gen = -1, beta_interv=-1):
    print('current sim ' + str(current))
    # With intervention into the simulation code
    sim = event_driven.Simulation(1000000, G, Lambda, Gamma, pos)
    sim.run_sim(intervention_gen, beta_interv)
    sm_matrix = sim.generate_matrix_gen(20)
    print('total timesteps', sim.total_num_timesteps)
    return sm_matrix


def generate_graph(N, deg_dist):
    number_of_nodes = N * np.array(deg_dist)
    degree_sequence = []
    for i in range(int(math.floor(len(number_of_nodes)))):
        number_with_that_degree = number_of_nodes[i]
        for k in range(int(math.floor(number_with_that_degree))):
            degree_sequence.append(i)
    # z = [5, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1]
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
    pos = nx.spring_layout(G)
    # nx.draw_networkx_nodes(G, pos=pos, with_labels=True)
    # nx.draw_networkx_labels(G, pos=pos, with_labels=True)
    # nx.draw_networkx_edges(G, pos=pos)
    # plt.show()
    # Check:
    inferred_degree_dist = np.array(nx.degree_histogram(G)) / N
    # print('Inferred equals given degree distribution: ', inferred_degree_dist == deg_dist)
    return G, pos
