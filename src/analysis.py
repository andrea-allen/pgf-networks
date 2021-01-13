import networkx as nx
import numpy as np
import math
from src import event_driven
import matplotlib.pyplot as plt
import time

def graph_infection_size_distribution_by_gen(list_of_gens, x_lim, filepath, filename, intervention_filepath=None, intervention_filename=None):
    intervention_comparison_true = False
    if intervention_filepath is not None:
        intervention_comparison_true = True
    color_key = {}
    colors = ['blue', 'red', 'orange', 'black', 'pink', 'teal', 'magenta', 'purple', 'green',
                 'aqua', 'brown', 'gold']
    for i in range(len(list_of_gens)):
        gen = list_of_gens[i]
        color = np.random.choice(colors)
        color_key[gen] = color

    data_no_intervention = np.loadtxt(filepath+filename, delimiter=',')

    for gen in list_of_gens:
        time_series = data_no_intervention[gen][2:x_lim]
        plt.plot(time_series, label='$g=$' + str(gen), color=color_key[gen], alpha=0.5, lw=1)

    if intervention_comparison_true:
        data_intervention = np.loadtxt(intervention_filepath+intervention_filename, delimiter=',')
        for gen in list_of_gens:
            time_series_int = data_intervention[gen][2:x_lim]
            plt.plot(time_series_int, color=color_key[gen], alpha=0.75, ls='--', lw=1)

    plt.legend(loc='upper right')
    plt.xlabel('number infected $s$ at generation $g$')
    plt.ylabel('$p_s^g$')
    # plt.title(
    #     'Effects of simulations with intervention from $T=.2$ to $T=.1$ at $g=3$ on Binomial Degree Distribution Network')
    plt.semilogy()
    plt.ylim(.0001, .1)
    # plt.title('Created from saved data')
    # plt.savefig('p_s_g_distribution_intervention.png')
    plt.show()

    # TODO add an inset plot option with total number of infections?

# def graph_phase_space():
    
