import networkx as nx
import numpy as np
import math
import event_driven
import matplotlib.pyplot as plt


def run():
    print('running')
    # simulate(0, 0, 0, 0)
    results = ensemble(200)
    for i in range(14):
        fig, ax = plt.subplots()
        ax.imshow(results[i], cmap='autumn')
        ax.invert_yaxis()
        plt.xlabel('s nodes infected')
        plt.ylabel('m nodes infected in gen g')
        # plt.legend(loc='upper right')
        plt.title('Generation g='+str(i))
        plt.savefig('gen'+str(i)+'samplefig.png')
        plt.show()


def ensemble(num_sims = 10):
    N = 30
    gen_results = np.zeros(((int(N/2), N,N)))
    for i in range(num_sims):
        sm_matrix = simulate(0, 0, 0, 0)
        for g in range(len(sm_matrix[0])):
            gen_s = int(sm_matrix[1][g])
            gen_m = int(sm_matrix[0][g])
            gen_results[g][gen_m][gen_s] += 1

    # averaging:
    for gen in range(int(N/2)):
        gen_matrix = gen_results[gen]
        gen_matrix = gen_matrix / num_sims
        gen_results[gen] = gen_matrix
    return gen_results
    #Want: One set of matrices per simulation:
# #   One matrix per gen g, matrix=sxm, s=number infected from 0, N, m=number infected during gen g: 0,N
    #Simulation_Data = []
    #for gen g in gens:
    # append(Matrix(number_tagged_in_g, number_tagged_in_g-number_tagged_in_g-1)
    # Can put the matrix data in the simulation class, then it spits out the whole matrix, per g? Or just have a self function on it
    # given all the matrices, take the sum of the set of matrices and divide by number of simulations?

def simulate(A, beta, gamma, num_times):
    degree_dist = [0, .20, .20, .15, .05, .35, .02, .02, .01]
    N = 30
    G, pos = generate_graph(N, degree_dist)
    Lambda = np.zeros((N, N))
    Gamma = np.zeros(N)
    for i in range(N):
        Gamma[i] = .05
        for j in range(N):
            Lambda[i][j] = .5
    sim = event_driven.Simulation(N*5, G, Lambda, Gamma, pos)
    sim.run_sim()
    sm_matrix = sim.generate_matrix_gen()
    return sm_matrix

def generate_graph(N, deg_dist):
    number_of_nodes = N*np.array(deg_dist)
    degree_sequence = []
    for i in range(int(math.floor(len(number_of_nodes)))):
        number_with_that_degree = number_of_nodes[i]
        for k in range(int(math.floor(number_with_that_degree))):
            degree_sequence.append(i)
    # z = [5, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1]
    graphical = nx.is_graphical(degree_sequence)
    # print('Degree sequence is graphical: ', graphical)
    if not graphical:
        print('Adding node of degree 1')
        degree_sequence.append(1)
    # print("Configuration model")
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
    # print('Inferred equals given degree distribution: ', inferred_degree_dist == deg_dist)
    return G, pos