import networkx as nx
import numpy as np
import math
import event_driven
import matplotlib.pyplot as plt


def run():
    print('running')
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

    results = ensemble(20)
    # TODO ^ fix this
    still_has_next = True
    while still_has_next:
        for i in range(len(results) - 1):
            fig, ax = plt.subplots()
            ax.imshow(results[i][:100][:, :100], cmap='autumn')
            ax.invert_yaxis()
            plt.xlabel('s nodes infected')
            plt.ylabel('m nodes infected in gen g')
            # plt.legend(loc='upper right')
            plt.title('Generation g=' + str(i))
            plt.savefig('gen' + str(i) + 'samplefig_noel.png')
            plt.show()
            still_has_next = np.sum(results[i + 1]) > 0


def power_law_degree_distrb():
    degree_dist = np.zeros(40)
    for k in range(1, len(degree_dist)):
        p_k = (k ** (-2)) * (math.e ** (-k / 5))
        degree_dist[k] = p_k
    return degree_dist


def ensemble(num_sims=10, N=1000):
    # ensemble of simulations to produce s,m phase space
    gen_results = np.zeros((int(N / 2), N, N))
    for i in range(num_sims):
        sm_matrix = simulate_noel(N, 0, 0, 0, 0)
        for g in range(len(sm_matrix[0])):
            gen_s = int(sm_matrix[1][g])
            gen_m = int(sm_matrix[0][g])
            gen_results[g][gen_m][gen_s] += 1

    # averaging:
    for gen in range(int(N / 2)):
        gen_matrix = gen_results[gen]
        gen_matrix = gen_matrix / num_sims
        gen_results[gen] = gen_matrix
    return gen_results
    # Want: One set of matrices per simulation:


def outbreak_size_distrb_per_gen(degree_distrb, num_sims=10, N=1000):
    s_sizes = np.arange(N)
    outbreak_size_distrb_per_gen_matrix = np.zeros((100, N))
    G, pos = generate_graph(N, degree_distrb)
    N = len(G.nodes())
    Lambda = np.zeros((N, N))
    Gamma = np.zeros(N)
    for n in range(N):
        Gamma[n] = .001
        for j in range(N):
            Lambda[n][j] = .8
    for i in range(num_sims):
        if i % 1000 == 0:
            G, pos = generate_graph(N, degree_distrb)
            N_resized = len(G.nodes())
            Lambda = np.zeros((N_resized, N_resized))
            Gamma = np.zeros(N_resized)
            for n in range(N_resized):
                Gamma[n] = .001
                for j in range(N_resized):
                    Lambda[n][j] = .8
        sm_matrix = simulate_noel(G, pos, Lambda, Gamma, i)
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
    sm_matrix = sim.generate_matrix_gen()
    return sm_matrix


def simulate_noel(G, pos, Lambda, Gamma, current):
    print('current sim ' + str(current))
    sim = event_driven.Simulation(1000000, G, Lambda, Gamma, pos)
    sim.run_sim()
    sm_matrix = sim.generate_matrix_gen()
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
        # print('Adding node of degree 1')
        degree_sequence.append(1)
    # print("Configuration model")
    G = nx.configuration_model(degree_sequence)
    # Remove self-loops and parallel edges
    # DONT DO self loops for now, its slowing evertyhing down, maybe audit later
    # try:
    #     G.remove_edges_from(nx.selfloop_edges(G))
    # except RuntimeError:
    #     print('No self loops to remove')
    pos = nx.spring_layout(G)
    # nx.draw_networkx_nodes(G, pos=pos, with_labels=True)
    # nx.draw_networkx_labels(G, pos=pos, with_labels=True)
    # nx.draw_networkx_edges(G, pos=pos)
    # plt.show()
    # Check:
    inferred_degree_dist = np.array(nx.degree_histogram(G)) / N
    # print('Inferred equals given degree distribution: ', inferred_degree_dist == deg_dist)
    return G, pos
