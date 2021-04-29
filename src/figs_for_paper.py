import time

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy import stats
import networkx as nx

from analysis import ensemble
from src import degree_distributions
from src import pgf_formalism
from src import plotting_util

import epintervene


def erdos_renyi_exp(file_root='poiss_T8_10k_q_1_gamma1_g_over_b',num_sims=50, num_nodes=10000, kill_by=None):
    ## Make erdos renyi network
    k_mean_degree = 2.5
    # k_mean_degree = 5
    er_degree_dist = degree_distributions.binomial_degree_distb(400, k_mean_degree)


    ## set up paramters
    beta = 0.004
    gamma = 0.001
    T = 0.8


    ## quantities to track:
    ## dist of cmulative infections by gen, dist of cum infections by time (g/beta*q?), total and active gens, diff of waiting times
    # file_root = 'poiss_T8_1k_q_1_gamma1_g_over_b'
    print('Running erdos renyi experiment')
    s_time = time.time()
    ensemble.run_ensemble_with_time_distribution(er_degree_dist, f'../data/testing/{file_root}', num_sims=num_sims,
                                                 num_nodes=num_nodes, init_T=T, recover_rate=gamma, kill_by=kill_by)
    print('Done')
    print(f'time all was {time.time() - s_time}')


def power_law_exp(file_root='plaw_T8_10k_q_1_gamma1_g_over_b',num_sims=50, num_nodes=10000, active_gen_sizes_on=False):
    ## Make power law network
    pl_degree_dist = degree_distributions.power_law_degree_distrb(400)
    power_law_q2 = degree_distributions.power_law_degree_distrb(400, mu=10)


    ## set up paramters
    beta = 0.004
    gamma = 0.001
    T = 0.8


    ## quantities to track:
    ## dist of cmulative infections by gen, dist of cum infections by time (g/beta*q?), total and active gens, diff of waiting times
    # file_root = 'poiss_T8_1k_q_1_gamma1_g_over_b'
    print('Running power law experiment')
    ensemble.run_ensemble_with_time_distribution(power_law_q2, f'../data/testing/{file_root}', num_sims=num_sims,
                                                 num_nodes=num_nodes, init_T=T, recover_rate=gamma, active_gen_sizes_on=active_gen_sizes_on)
    print('Done')

def tree_exp(file_root='tree_T8_10k', num_sims=50, num_nodes=8191):
    ## Make tree network
    # TODO
    tree_graph = nx.generators.balanced_tree(2, 12)
    tree_dd = np.array(nx.degree_histogram(tree_graph))
    tree_dd = tree_dd/len(nx.nodes(tree_graph))


    ## set up paramters
    beta = 0.004
    gamma = 0.001
    T = 0.8


    ## quantities to track:
    ## dist of cmulative infections by gen, dist of cum infections by time (g/beta*q?), total and active gens, diff of waiting times
    # file_root = 'poiss_T8_1k_q_1_gamma1_g_over_b'
    print('Running power law experiment')
    ensemble.run_ensemble_with_time_distribution(tree_dd, f'../data/testing/{file_root}', num_sims=num_sims, num_nodes=num_nodes, init_T=T, recover_rate=gamma)
    print('Done')

def results_plots(file_root='poiss_T8_10k_q_1_gamma1_g_over_b', q_degree=1, active_gen_sizes_on=False):
    plt.figure('Generations')
    plt.title('Generational time distribution of cumulative infections')
    gen_emergence = np.loadtxt(f'../data/testing/{file_root}_gen_emergence_times.txt', delimiter=',')
    plotting_util.graph_infection_size_distribution_by_gen([2, 4, 6, 10, 20], 200, f'../data/testing/{file_root}_generational.txt', gen_emergence_times=gen_emergence)
    # vlines doesn't really make sense here since the x axis is node quantity, but, can use as labels
    # plt.vlines(gen_emergence, ymin=0, ymax=1, colors='orange', alpha=0.5, linestyles=':')
    plt.tight_layout()
    plt.savefig(f'../data/{file_root}_gens.png')
    plt.figure('Time')
    t_buckets = np.loadtxt(f'../data/testing/{file_root}_time_distribution_values.txt', delimiter=',')
    t_distribution = np.loadtxt(f'../data/testing/{file_root}_time_distribution.txt', delimiter=',')
    colors = ['red', 'orange', 'green', 'blue', 'goldenrod', 'purple', 'pink', 'teal', 'black', 'gold', 'chocolate',
              'dodgerblue', 'darkslategray', 'mediumorchid', 'magenta', 'tomato', 'midnightblue',
              'cadetblue', 'crimson']
    standrd_colors = ['black', 'blue', 'red', 'purple', 'orange', 'brown']

    for x in [3, 5, 7, 11, 21]:
        plt.plot(t_distribution[x][2:200], label=f't={int(t_buckets[x])}, Expected gen (txBeta) = {int(x)}', color=standrd_colors[0], alpha=0.7)
        colors.remove(colors[0])
        standrd_colors.remove(standrd_colors[0])
    plt.legend(loc='upper right')
    plt.semilogy()
    plt.xlabel('Cumulative number nodes infected')
    plt.ylabel('Probability')
    plt.title('Clock time distribution of cumulative infections')
    plt.tight_layout()
    plt.savefig(f'../data/{file_root}_time.png')
    plt.figure("Generational emergence vs expectation")
    plt.title('Typical clock time of generation emergence vs expected clock time')
    max_len = min(len(t_buckets), len(gen_emergence))
    # max_len = 60
    max_len = 21
    plt.scatter(t_buckets[:max_len], gen_emergence[:max_len])
    plt.xlabel('Expected time of emergence')
    plt.ylabel('Actual time of emergence')
    lin_reg_result = stats.linregress(t_buckets[:max_len], gen_emergence[:max_len])
    x_vals = np.arange(0, 15000)
    x_vals = np.arange(0, 2000)
    y_fit = lin_reg_result.intercept + lin_reg_result.slope * x_vals
    plt.plot(x_vals, y_fit, color='red', label=f'Linear fit of $y={np.round(lin_reg_result.slope, 2)}$')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(f'../data/{file_root}_time_fit.png')
    plt.show()

    expected_time = np.zeros(100)
    beta = .004
    q_gcc = 1.5
    q_main = q_degree
    for g in range(1, 100):
        expected_time[g] = g/(q_main*beta)
    plt.scatter(expected_time[:20], gen_emergence[:20], label='q_main')
    plt.show()
    plt.scatter(np.arange(40), expected_time[:40], label='expected, g/(q*beta) increments')
    plt.scatter(np.arange(40), gen_emergence[:40], label='actual average time of emergence')
    plt.title(f'Modified power law with <k>=1.5, \n <q>={q_main}')
    plt.legend(loc='upper left')
    plt.show()
    # expected_time = np.zeros(100)
    # for g in range(1, 100):
    #     expected_time[g] = np.sum(expected_time[:g]) + 1/((q_gcc*beta))
    # plt.scatter(expected_time, gen_emergence, label='q_gcc')
    # plt.legend(loc='upper left')
    # plt.show()

    # correlation = plotting_util.gen_vs_time_correlation(np.loadtxt(f'../data/testing/{file_root}_generational.txt', delimiter=','), t_distribution)
    # plt.scatter(np.arange(len(t_distribution))[2:], correlation[2:])
    # plt.show()
    fig, ax1 = plt.subplots(figsize=(14, 7))

    active_gens_ts = np.loadtxt(f'../data/testing/{file_root}_active_gen_ts.txt', delimiter=',')
    total_gens_ts = np.loadtxt(f'../data/testing/{file_root}_total_gen_ts.txt', delimiter=',')
    timeseries_vals = np.loadtxt(f'../data/testing/{file_root}_ts_vals_normalized.txt', delimiter=',')
    ax1.plot(timeseries_vals, active_gens_ts, label='active generations')
    ax1.plot(timeseries_vals, total_gens_ts, label='total generations')
    # ax1.set_xticks(
    #     [timeseries_vals[0], timeseries_vals[50], timeseries_vals[100], timeseries_vals[150], timeseries_vals[200],
    #      timeseries_vals[250]])
    q_main = q_degree
    start_color = (.1, .2, .5)
    for g in range(1, 15):
        expected_time[g] = g / (q_main * beta)
    for g in [2, 6, 10, 12, 18, 30]:
        ax1.vlines(expected_time[g], ymin=0, ymax=max(total_gens_ts), color=(.1, .2, .3 + g / 50), ls='--')
        ax1.vlines(gen_emergence[g], ymin=0, ymax=max(total_gens_ts), color=(.1, .2, .3 + g / 50))
    ax1.legend(loc='upper left')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Number active generations at time t')
    # plt.show()
    if not active_gen_sizes_on:
        plt.show()

    if active_gen_sizes_on:
        right, bottom, width, height = [0.135, 0.15, 0.75, 0.2]
        ax2 = fig.add_axes([right, bottom, width, height])
        # ax2.set_xticks([timeseries_vals[0], timeseries_vals[50], timeseries_vals[100], timeseries_vals[150], timeseries_vals[200], timeseries_vals[250]])
        # TODO make it an actual inset
        average_active_sizes = np.loadtxt(f'../data/testing/{file_root}_active_gen_sizes_ts.txt', delimiter=',')
        active_gens_inset(timeseries_vals, average_active_sizes, ax2)

    # Same plot for gens and time
    plt.figure('Generations')
    plt.title('Generational time distribution of cumulative infections')
    gen_emergence = np.loadtxt(f'../data/testing/{file_root}_gen_emergence_times.txt', delimiter=',')
    plt.subplot(2,1,1)
    plotting_util.graph_infection_size_distribution_by_gen([4, 6, 10, 20], 200,
                                                           f'../data/testing/{file_root}_generational.txt',
                                                           gen_emergence_times=gen_emergence)
    # vlines doesn't really make sense here since the x axis is node quantity, but, can use as labels
    # plt.vlines(gen_emergence, ymin=0, ymax=1, colors='orange', alpha=0.5, linestyles=':')
    plt.xlabel('Cumulative infections')
    plt.ylabel('Probability')
    plt.ylim([.0001, .05])
    plt.tight_layout()
    t_buckets = np.loadtxt(f'../data/testing/{file_root}_time_distribution_values.txt', delimiter=',')
    t_distribution = np.loadtxt(f'../data/testing/{file_root}_time_distribution.txt', delimiter=',')
    colors = ['red', 'orange', 'green', 'blue', 'goldenrod', 'purple', 'pink', 'teal', 'black', 'gold', 'chocolate',
              'dodgerblue', 'darkslategray', 'mediumorchid', 'magenta', 'tomato', 'midnightblue',
              'cadetblue', 'crimson']
    standrd_colors = ['black', 'blue', 'red', 'purple', 'orange', 'brown']

    # TODO fix axes, labels, etc. Then write captions tonight 530-7
    plt.subplot(2,1,2)
    for x in [5, 7, 11, 21]:
        if x==5: #first entry
            plt.plot(t_distribution[x][2:200], label=f'Cumulative infections by t={int(t_buckets[x])}={x}/$q\\beta$',
                     color=standrd_colors[0], alpha=0.7, ls='-', lw=1)
        else:
            plt.plot(t_distribution[x][2:200], label=f't={int(t_buckets[x])}={x}/$q\\beta$',
                     color=standrd_colors[0], alpha=0.7, ls='-', lw=1)
        # colors.remove(colors[0])
        standrd_colors.remove(standrd_colors[0])
    plt.semilogy()
    plt.legend(loc='upper right')
    plt.xlabel('Cumulative infections')
    plt.ylim([.0001, .05])
    plt.ylabel('Probability')
    plt.tight_layout()
    plt.show()

    ## Trying to make figure 3 again
    if active_gen_sizes_on:
        gens_to_display_lines = [2, 6, 10, 14, 18, 22, 26, 30]
        gens_to_display_curves = [2, 4, 6, 8, 10, 12]
        ax2 = plt.subplot(212)
        # ax1 = plt.subplot(211, sharex=ax2)
        ax1 = plt.subplot(211)
        max_x_val = int(.8*len(timeseries_vals))

        active_gens_ts = np.loadtxt(f'../data/testing/{file_root}_active_gen_ts.txt', delimiter=',')
        total_gens_ts = np.loadtxt(f'../data/testing/{file_root}_total_gen_ts.txt', delimiter=',')
        timeseries_vals = np.loadtxt(f'../data/testing/{file_root}_ts_vals_normalized.txt', delimiter=',')
        # plt.subplot(2, 1, 1, sharex=True)
        ax1.plot(timeseries_vals[:max_x_val], active_gens_ts[:max_x_val], label='active \n generations', color='blue', ls='--')
        ax1.plot(timeseries_vals[:max_x_val], total_gens_ts[:max_x_val], label='total \n generations', color='blue', ls='-')
        # ax1.set_xticks(
        #     [timeseries_vals[0], timeseries_vals[50], timeseries_vals[100], timeseries_vals[150], timeseries_vals[200],
        #      timeseries_vals[250]])
        q_main = q_degree
        start_color = (.1, .2, .5)
        for g in range(1, 15):
            expected_time[g] = g / (q_main * beta)
        for g in gens_to_display_lines:
            if g == 2: #just do a label for one pair
                ax1.vlines(expected_time[g], ymin=0, ymax=max(total_gens_ts), color=(.1, .2, .3 + g / 50), ls=':')
                           # ,alpha=0.5, label='Expected time of \n emergence $\\frac{g}{q\\beta}$')
                ax1.vlines(gen_emergence[g], ymin=0, ymax=max(total_gens_ts), color=(.1, .2, .3 + g / 50), alpha=0.5)
                           # ,label='Average empirical time \n of emergence')
            else:
                ax1.vlines(expected_time[g], ymin=0, ymax=max(total_gens_ts), color=(.1, .2, .3 + g / 50), ls=':', alpha=0.5)
                ax1.vlines(gen_emergence[g], ymin=0, ymax=max(total_gens_ts), color=(.1, .2, .3 + g / 50), alpha=0.5)
        ax1.legend(loc='upper right')
        ax1.set_xticks(list([expected_time[g] for g in gens_to_display_lines]))
        x_tick_labels = list([f'$\\frac{{{g}}}{{q\\beta}}$' for g in gens_to_display_lines])
        x_tick_labels[0] = f'$t = \\frac{{g}}{{q\\beta}}$'
        ax1.set_xticklabels(x_tick_labels)
        ax1.xaxis.tick_top()
        ax1.xaxis.set_label_position('top')
        ax1.tick_params(axis='x', labelrotation=0, labelsize=12)
        # ax1.set_xticks([])
        # ax1.set_xlabel('Time')
        ax1.set_ylabel('Number active generations')
        # plt.show()
        # plt.subplot(2, 1, 2, sharex=True)
        # ax2.set_xticks([timeseries_vals[0], timeseries_vals[50], timeseries_vals[100], timeseries_vals[150], timeseries_vals[200], timeseries_vals[250]])
        # TODO make it an actual inset
        average_active_sizes = np.loadtxt(f'../data/testing/{file_root}_active_gen_sizes_ts.txt', delimiter=',')
        for g in gens_to_display_curves:
            ax2.plot(timeseries_vals[:max_x_val], average_active_sizes.T[g][:max_x_val], color=(.1, .2, .3 + g / 20))
            x_position_val = np.argmax(average_active_sizes.T[g])
            max_y_val = average_active_sizes.T[g][x_position_val]
            ax2.text(timeseries_vals[x_position_val]+200, max_y_val-1, f'g {g}', horizontalalignment='left')
        for g in gens_to_display_lines:
            if g == 2: #just do a label for one pair
                ax2.vlines(expected_time[g], ymin=0, ymax=np.max(average_active_sizes), color=(.1, .2, .3 + g / 50), ls=':',
                           alpha=0.5, label='Expected time \n of emergence $\\frac{g}{q\\beta}$')
                ax2.vlines(gen_emergence[g], ymin=0, ymax=np.max(average_active_sizes), color=(.1, .2, .3 + g / 50), alpha=0.5,
                           label='Average empirical \n time of gen $g$ \n emergence')
            else:
                ax2.vlines(expected_time[g], ymin=0, ymax=np.max(average_active_sizes), color=(.1, .2, .3 + g / 50), alpha=0.5, ls=':')
                ax2.vlines(gen_emergence[g], ymin=0, ymax=np.max(average_active_sizes), color=(.1, .2, .3 + g / 50), alpha=0.5)
        ax2.legend(loc='upper right')
        # plt.xlabel('time')
        ax2.set_xlabel('Time')
        ax2.set_xticks(np.arange(0, 4001, 500))
        x_tick_labels_time = list(np.arange(0, 4001, 500))
        x_tick_labels_time[0] = '$t=0$'
        ax2.set_xticklabels(x_tick_labels_time)
        ax2.tick_params(axis='x', labelrotation=-20)
        ax2.set_ylabel('Number active nodes')
        plt.show()

def active_gens_inset(timeseries, average_active_sizes, ax2):
    #TODO active gens
    for g in range(1, 12, 3):
        ax2.plot(timeseries, average_active_sizes.T[g], label=f'generation {g}')
    ax2.legend(loc='upper right')
    ax2.set_xlabel('time')
    ax2.set_ylabel('average number active nodes')
    plt.show()

def combine_data():
    results_1 = np.loadtxt(f'../data/testing/erdos_renyi_10k_10ksims_pt1_generational.txt', delimiter=',')
    results_2 = np.loadtxt(f'../data/testing/erdos_renyi_10k_10ksims_pt2_generational.txt', delimiter=',')
    results_3 = np.loadtxt(f'../data/testing/erdos_renyi_10k_10ksims_pt3_generational.txt', delimiter=',')
    results_4 = np.loadtxt(f'../data/testing/erdos_renyi_10k_10ksims_pt4_generational.txt', delimiter=',')
    results_5 = np.loadtxt(f'../data/testing/erdos_renyi_10k_10ksims_pt5_generational.txt', delimiter=',')
    results_6 = np.loadtxt(f'../data/testing/erdos_renyi_10k_10ksims_pt6_generational.txt', delimiter=',')
    combo_results = np.array(results_1) + np.array(results_2) + np.array(results_3) + np.array(results_4) + np.array(results_5) + np.array(results_6)
    combo_results = combo_results/6
    # np.savetxt(f'../data/testing/erdos_renyi_10k_60ksims_combo_generational.txt', combo_results, delimiter=',')
    print(combo_results[1][:10])


