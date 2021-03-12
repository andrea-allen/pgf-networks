import numpy as np
import matplotlib.pyplot as plt


def graph_infection_size_distribution_by_gen(list_of_gens, x_lim, filepath, filename, intervention_filepath=None,
                                             intervention_filename=None):
    intervention_comparison_true = False
    if intervention_filepath is not None:
        intervention_comparison_true = True
    color_key = {}
    colors = ['red', 'orange', 'green', 'blue', 'purple', 'teal', 'black', 'gold', 'chocolate',
              'dodgerblue', 'darkslategray', 'mediumorchid', 'magenta', 'tomato', 'midnightblue',
              'cadetblue', 'crimson']
    for i in range(len(list_of_gens)):
        gen = list_of_gens[i]
        color = np.random.choice(colors)
        colors.remove(color)
        color_key[gen] = color

    data_no_intervention = np.loadtxt(filepath + filename, delimiter=',')

    for gen in list_of_gens:
        time_series = data_no_intervention[gen][2:x_lim]
        plt.plot(time_series, label='$g=$' + str(gen), color=color_key[gen], alpha=0.5, lw=1)

    if intervention_comparison_true:
        data_intervention = np.loadtxt(intervention_filepath + intervention_filename, delimiter=',')
        for gen in list_of_gens:
            time_series_int = data_intervention[gen][2:x_lim]
            plt.plot(time_series_int, color=color_key[gen], alpha=0.75, ls='--', lw=1)

    plt.legend(loc='upper right')
    plt.xlabel('number infected $s$ at generation $g$')
    plt.ylabel('$p_s^g$')
    plt.semilogy()
    plt.ylim(.0001, .1)
