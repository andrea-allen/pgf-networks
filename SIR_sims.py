import networkx as nx
import numpy as np
import math
import event_driven


def run():
    print('running')
    simulate(0, 0, 0, 0)

def simulate(A, beta, gamma, num_times):
    degree_dist = [0, .20, .20, .15, .05, .35, .02, .02, .01]
    N = 40
    G, pos = generate_graph(N, degree_dist)
    A = nx.adjacency_matrix(G)

    Lambda = np.zeros((N, N))
    Gamma = np.zeros(N)
    for i in range(N):
        Gamma[i] = .05
        for j in range(N):
            Lambda[i][j] = .5
    sim = event_driven.Simulation(100, G, Lambda, Gamma, pos)
    sim.run_sim()
    print(len(sim.has_been_infected_labels))


def generate_graph(N, deg_dist):
    number_of_nodes = N*np.array(deg_dist)
    degree_sequence = []
    for i in range(int(math.floor(len(number_of_nodes)))):
        number_with_that_degree = number_of_nodes[i]
        for k in range(int(math.floor(number_with_that_degree))):
            degree_sequence.append(i)
    # z = [5, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1]
    graphical = nx.is_graphical(degree_sequence)
    print('Degree sequence is graphical: ', graphical)
    if not graphical:
        print('Adding node of degree 1')
        degree_sequence.append(1)
    print("Configuration model")
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
    inferred_degree_dist = np.array(nx.degree_histogram(G))/N
    print('Inferred equals given degree distribution: ', inferred_degree_dist == deg_dist)
    return G, pos