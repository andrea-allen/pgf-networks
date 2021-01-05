import networkx as nx
import numpy as np
import math
from src import event_driven
from src import SIR_sims
import matplotlib.pyplot as plt
import time

def run_sandbox():
    s_sizes, data = run_single_simulation()
    color_key = {2: 'blue', 6: 'red', 11: 'orange', 18:'black'}
    for gen in [2, 6, 11]:
        time_series = data[gen][2:200]
        # time_series_int = data_int[gen][2:200]
        plt.plot(time_series, label='$g=$' + str(gen), color=color_key[gen], alpha=0.5, lw=1)
        # plt.plot(time_series_int, color=color_key[gen], alpha=0.75, ls='--', lw=1)
    plt.legend(loc='upper right')
    plt.xlabel('number infected $s$ at generation $g$')
    plt.ylabel('$p_s^g$')
    # plt.title('Effects of simulations with intervention from $T=.2$ to $T=.1$ at $g=3$ on Binomial Degree Distribution Network')
    plt.semilogy()
    plt.ylim(.0001, .1)
    # plt.title('Created from saved data')
    # plt.savefig('p_s_g_distribution_intervention.png')
    plt.show()

def run_single_simulation():
    degree_distrb = SIR_sims.power_law_degree_distrb()
    gamma = 0.1
    T = 0.8
    N = 10000
    beta = -(gamma * T) / (T - 1)
    s_sizes = np.arange(N)
    outbreak_size_distrb_per_gen_matrix = np.zeros((100, N))
    G, pos = SIR_sims.generate_graph(N, degree_distrb)
    N = len(G.nodes())
    Lambda = np.zeros((N, N))
    Gamma = np.zeros(N)
    # Construct recovery values and transmissibility matrix
    for n in range(N):
        Gamma[n] = gamma
        for j in range(N):
            Lambda[n][j] = beta
    # Run simulations, new graph every 500 simulations
    start_time = time.time()
    results = SIR_sims.simulate(G, pos, Lambda, Gamma, 1)  # results for time series of s and m nodes infected
    print("--- %s seconds to run simulation---" % (time.time() - start_time))
    for g in range(len(results[0])):
        gen_s = int(results[1][g])  # Second vector is total infections over time in s
        try:
            outbreak_size_distrb_per_gen_matrix[g][gen_s] += 1  # Add to mass of that result
        except IndexError:
            print('Index error for g: ', g, ', gen_s: ', gen_s)
            continue
    return s_sizes, outbreak_size_distrb_per_gen_matrix