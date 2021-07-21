import time

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy import stats

from analysis import ensemble
from src import degree_distributions
from src import plotting_util
import seaborn as sns


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
                                                 num_nodes=num_nodes, init_T=T, recover_rate=gamma,
                                                 active_gen_sizes_on=active_gen_sizes_on, kill_by=16)
    print('Done')

def random_vax_exp(file_root, num_sims=50, num_nodes=10000):
    ## Make erdos renyi network
    k_mean_degree = 3.04
    # k_mean_degree = 5
    power_law_q2 = degree_distributions.power_law_degree_distrb(400, mu=10)


    ## set up paramters
    beta = 0.004
    gamma = 0.001
    T = 0.8


    ## quantities to track:
    ## dist of cmulative infections by gen, dist of cum infections by time (g/beta*q?), total and active gens, diff of waiting times
    # file_root = 'poiss_T8_1k_q_1_gamma1_g_over_b'
    print('Running power law intervention experiment')
    s_time = time.time()
    ensemble.run_ensemble_intervention_effects(power_law_q2, f'../data/testing/{file_root}', num_sims=num_sims,
                                               num_nodes=num_nodes, init_T=T, gen_intervene=4, T_intervene=0.0,
                                               recover_rate=gamma, prop_reduced=0.5, intervention_type="random",
                                               run_regular=False)
    print('Done')
    print(f'time all was {time.time() - s_time}')

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
    ## FIGURE FOR PAPER
    if active_gen_sizes_on:
        beta = .004
        time_emergence_plot(file_root, q_degree, beta, save=True)

    ## NOT FOR PAPER
    gen_emergence = np.loadtxt(f'../data/{file_root}_gen_emergence_times.txt', delimiter=',')
    # vlines doesn't really make sense here since the x axis is node quantity, but, can use as labels
    # plt.vlines(gen_emergence, ymin=0, ymax=1, colors='orange', alpha=0.5, linestyles=':')
    t_buckets = np.loadtxt(f'../data/{file_root}_time_distribution_values.txt', delimiter=',')
    t_distribution = np.loadtxt(f'../data/{file_root}_time_distribution.txt', delimiter=',')
    colors = ['red', 'orange', 'green', 'blue', 'goldenrod', 'purple', 'pink', 'teal', 'black', 'gold', 'chocolate',
              'dodgerblue', 'darkslategray', 'mediumorchid', 'magenta', 'tomato', 'midnightblue',
              'cadetblue', 'crimson']
    standrd_colors = ['black', 'blue', 'red', 'purple', 'orange', 'brown']

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

    # Same plot for gens and time
    ## POSSIBLE FIGURE FOR PAPER
    plt.figure('Generations')
    plt.title('Generational time distribution of cumulative infections')
    gen_emergence = np.loadtxt(f'../data/{file_root}_gen_emergence_times.txt', delimiter=',')
    plt.subplot(2,1,1)
    plotting_util.graph_infection_size_distribution_by_gen([4, 6, 10, 20], 200,
                                                           f'../data/{file_root}_generational.txt',
                                                           gen_emergence_times=gen_emergence)
    # vlines doesn't really make sense here since the x axis is node quantity, but, can use as labels
    # plt.vlines(gen_emergence, ymin=0, ymax=1, colors='orange', alpha=0.5, linestyles=':')
    plt.xlabel('Cumulative infections')
    plt.ylabel('Probability')
    plt.ylim([.0001, .05])
    plt.tight_layout()
    t_buckets = np.loadtxt(f'../data/{file_root}_time_distribution_values.txt', delimiter=',')
    t_distribution = np.loadtxt(f'../data/{file_root}_time_distribution.txt', delimiter=',')
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





def time_emergence_plot(file_root, q_degree, beta, save=True):
    #TODO: bone colors?
    gens_to_display_lines = [2, 4, 6, 8, 10, 12, 14, 18, 20, 22]
    gens_to_display_lines = [2, 6, 10,  14, 18, 22, 26, 28]
    gens_to_display_curves = [4, 8, 10, 14]
    cmap = plt.get_cmap('bone')
    indices = np.linspace(0, cmap.N, int(max(gens_to_display_lines)*1.5) + 2)
    my_colors = [cmap(int(i)) for i in indices]
    standrd_colors = my_colors[:max(gens_to_display_lines)+2][::-1] #Shifted to avoid yellow
    ax2 = plt.subplot(212)
    # ax1 = plt.subplot(211, sharex=ax2)
    ax1 = plt.subplot(211)

    expected_time = np.zeros(100)
    for g in range(1, 100):
        expected_time[g] = g/(q_degree*beta)

    gen_emergence = np.loadtxt(f'../data/{file_root}_gen_emergence_times.txt', delimiter=',')
    active_gens_ts = np.loadtxt(f'../data/{file_root}_active_gen_ts.txt', delimiter=',')
    total_gens_ts = np.loadtxt(f'../data/{file_root}_total_gen_ts.txt', delimiter=',')
    timeseries_vals = np.loadtxt(f'../data/{file_root}_ts_vals_normalized.txt', delimiter=',')
    max_x_val = int(.35 * len(timeseries_vals))
    # plt.subplot(2, 1, 1, sharex=True)
    offset_for_top_plot = int(max_x_val/10)
    ax1.plot(timeseries_vals[:max_x_val - offset_for_top_plot], active_gens_ts[:max_x_val- offset_for_top_plot], label='active \n generations', color=standrd_colors[-1],
             ls='--')
    ax1.plot(timeseries_vals[:max_x_val - offset_for_top_plot], total_gens_ts[:max_x_val - offset_for_top_plot], label='total \n generations', color=standrd_colors[-1], ls='-')
    ax1.text(timeseries_vals[max_x_val - offset_for_top_plot]+5, active_gens_ts[max_x_val - offset_for_top_plot], 'active', horizontalalignment='left', fontsize=12)
    ax1.text(timeseries_vals[max_x_val - offset_for_top_plot]+5, total_gens_ts[max_x_val - offset_for_top_plot], 'total', horizontalalignment='left', fontsize=12)
    # ax1.set_xticks(
    #     [timeseries_vals[0], timeseries_vals[50], timeseries_vals[100], timeseries_vals[150], timeseries_vals[200],
    #      timeseries_vals[250]])
    ax1.set_xlim([0, timeseries_vals[max_x_val+1]])
    q_main = q_degree
    start_color = (.1, .2, .5)
    for g in range(1, 15):
        expected_time[g] = g / (q_main * beta)
    for i in range(len(gens_to_display_lines)):
        g = gens_to_display_lines[i]
        if g == 2:  # just do a label for one pair
            ax1.vlines(expected_time[g], ymin=0, ymax=max(total_gens_ts), color=standrd_colors[g], ls=':')
            # ,alpha=0.5, label='Expected time of \n emergence $\\frac{g}{q\\beta}$')
            ax1.vlines(gen_emergence[g], ymin=0, ymax=max(total_gens_ts), color=standrd_colors[g], alpha=0.5)
            # ,label='Average empirical time \n of emergence')
            expctd_time = expected_time[g]
            actl_time = gen_emergence[g]
            if expctd_time < actl_time:
                ax1.axvspan(expctd_time, actl_time, alpha=0.5, color='green')
            elif expctd_time > actl_time:
                ax1.axvspan(expctd_time, actl_time, alpha=0.5, color='yellow')
        else:
            ax1.vlines(expected_time[g], ymin=0, ymax=max(total_gens_ts), color=standrd_colors[g], ls=':',
                       alpha=0.5)
            ax1.vlines(gen_emergence[g], ymin=0, ymax=max(total_gens_ts), color=standrd_colors[g], alpha=0.5)
            expctd_time = expected_time[g]
            actl_time = gen_emergence[g]
            if expctd_time < actl_time:
                ax1.axvspan(expctd_time, actl_time, alpha=0.5, color='green')
            elif expctd_time > actl_time:
                ax1.axvspan(expctd_time, actl_time, alpha=0.5, color='yellow')
    ax1.text(expected_time[gens_to_display_lines[-1]]+4, 10, '$\\leftarrow$ expected time of gen', horizontalalignment='left', rotation=0, fontsize=12)
    ax1.text(gen_emergence[gens_to_display_lines[-1]]+4, 3, '$\\leftarrow$ empirical time of gen', horizontalalignment='left', rotation=0, fontsize=12)
    # ax1.legend(loc='upper right', frameon=False)
    ax1.set_xticks(list([expected_time[g] for g in gens_to_display_lines]))
    # x_tick_labels = list([f'$\\frac{{{g}}}{{q\\beta}}$' for g in gens_to_display_lines])
    # x_tick_labels[0] = f'$t = \\frac{{g}}{{q\\beta}}$'
    # ax1.set_xticklabels(x_tick_labels)
    ax1.xaxis.tick_top()
    ax1.xaxis.set_label_position('top')
    ax1.tick_params(axis='x', labelrotation=-45) #, labelsize=12)
    # ax1.set_xticks([])
    # ax1.set_xlabel('Time')
    ax1.set_ylabel('Number generations', fontsize=12)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(True)
    # plt.show()
    # plt.subplot(2, 1, 2, sharex=True)
    # ax2.set_xticks([timeseries_vals[0], timeseries_vals[50], timeseries_vals[100], timeseries_vals[150], timeseries_vals[200], timeseries_vals[250]])
    # TODO make it an actual inset
    average_active_sizes = np.loadtxt(f'../data/{file_root}_active_gen_sizes_ts.txt', delimiter=',')
    for i in range(len(gens_to_display_curves)):
        g = gens_to_display_curves[i]
        ax2.plot(timeseries_vals[:max_x_val - offset_for_top_plot], average_active_sizes.T[g][:max_x_val-offset_for_top_plot], color=standrd_colors[g])
        x_position_val = np.argmax(average_active_sizes.T[g]) #for specific positioning of labels
        static_position = max_x_val - int(max_x_val / 10) # for static positioning of curve labels
        max_y_val = average_active_sizes.T[g][static_position]
        ax2.text(timeseries_vals[static_position]+5, max_y_val, f'$g={{{g}}}$', horizontalalignment='left', fontsize=12)
    ax2.text(0, np.max(average_active_sizes)+80, 'gen $g$:', horizontalalignment='center', fontsize=12)
    for i in range(len(gens_to_display_lines)):
        g = gens_to_display_lines[i]
        if g == 2:  # just do a label for one pair
            ax2.vlines(expected_time[g], ymin=0, ymax=np.max(average_active_sizes), color=standrd_colors[g], ls=':',
                       alpha=0.5, label='Expected time \n of emergence $\\frac{g}{q\\beta}$')
            ax2.vlines(gen_emergence[g], ymin=0, ymax=np.max(average_active_sizes), color=standrd_colors[g],
                       alpha=0.5,
                       label='Average empirical \n time of gen $g$ \n emergence')
            expctd_time = expected_time[g]
            actl_time = gen_emergence[g]

            if expctd_time < actl_time:
                ax2.axvspan(expctd_time, actl_time, alpha=0.5, color='green', ymax=np.max(average_active_sizes))
                midpnt = expctd_time + (actl_time - expctd_time)/2
                ax2.text(midpnt, np.max(average_active_sizes) + 80, f'{g}', horizontalalignment='center',
                         rotation=0, fontsize=12)
            elif expctd_time > actl_time:
                ax2.axvspan(expctd_time, actl_time, alpha=0.5, color='yellow', ymax=np.max(average_active_sizes))
                midpnt = actl_time + (expctd_time - actl_time)/2
                ax2.text(midpnt, np.max(average_active_sizes) + 80, f'{g}', horizontalalignment='center',
                         rotation=0, fontsize=12)
        else:
            ax2.vlines(expected_time[g], ymin=0, ymax=np.max(average_active_sizes), color=standrd_colors[g],
                       alpha=0.5, ls=':')
            ax2.vlines(gen_emergence[g], ymin=0, ymax=np.max(average_active_sizes), color=standrd_colors[g],
                       alpha=0.5)
            # ax2.plot(gen_emergence[g], np.max(average_active_sizes)+120, 'o', color='black')
            expctd_time = expected_time[g]
            actl_time = gen_emergence[g]
            # Text for which generation the bar is for:
            # ax2.text(gen_emergence[g], np.max(average_active_sizes)+80, f'{g}', horizontalalignment='center', rotation=0)

            if expctd_time < actl_time:
                ax2.axvspan(expctd_time, actl_time, alpha=0.5, color='green', ymax=np.max(average_active_sizes))
                midpnt = expctd_time + (actl_time - expctd_time)/2
                ax2.text(midpnt, np.max(average_active_sizes) + 80, f'{g}', horizontalalignment='center',
                         rotation=0, fontsize=12)
            elif expctd_time > actl_time:
                ax2.axvspan(expctd_time, actl_time, alpha=0.5, color='yellow', ymax=np.max(average_active_sizes))
                midpnt = actl_time + (expctd_time - actl_time)/2
                ax2.text(midpnt, np.max(average_active_sizes) + 80, f'{g}', horizontalalignment='center',
                         rotation=0, fontsize=12)
    # ax2.legend(loc='upper right', frameon=False)
    # plt.xlabel('time')
    ax2.set_xlabel('Time (unitless)', fontsize=12)
    # ax2.set_xticks(np.arange(0, int(timeseries_vals[:max_x_val][-1]), 500))
    ax2.set_xticks(list([gen_emergence[g] for g in gens_to_display_lines]))
    # x_tick_labels_time = list(np.arange(0, int(timeseries_vals[:max_x_val][-1]), 500))
    # x_tick_labels_time[0] = '$t=0$'
    # ax2.set_xticklabels(x_tick_labels_time)
    ax2.tick_params(axis='x', labelrotation=-45)
    ax2.set_ylabel('Number active nodes in gen', fontsize=12)
    ax2.set_xlim([0, timeseries_vals[max_x_val+1]])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(True)
    ax2.spines['left'].set_visible(True)
    plt.tight_layout()
    if save:
        plt.savefig('time_plot_color.png')
    plt.savefig('time_plot_svg.svg', format='svg')
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

def combine_plaw_results():
    results_june = np.loadtxt(f'../data/testing/plaw_06_03_generational.txt', delimiter=',') #50,000 simulations
    results_may = np.loadtxt(f'../data/testing/plaw_05_25_generational.txt', delimiter=',') #75000 simulations
    combo_results = (np.array(results_june) * 50) + (np.array(results_may) * 75)
    combo_results = combo_results / (50+75)
    np.savetxt(f'../data/paper/plaw_combo_125k_generational.txt', combo_results, delimiter=',')
    print(combo_results[0])

def network_drawing():
    # nodes = list(np.arange(8))
    edges = [(0,1), (0,2), (0,3), (1,4), (1,5), (3,6), (3,7), (5,8), (6,7), (0,9), (9,10), (9,11)]
    G = nx.Graph()
    # G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    fig, axs = plt.subplots(2, 2)
    pos1 = nx.spring_layout(G)

    node_labels = {0:'0', 1:'', 2:'', 3:'', 4:'', 5:'', 6:'', 7:'', 8:'', 9:'', 10:'', 11:''}
    cmap = sns.color_palette('colorblind', n_colors=4)
    susceptible_color = cmap[1]
    infected_color = cmap[2]
    node_colors_step1 = [infected_color, susceptible_color, susceptible_color,susceptible_color,susceptible_color
        ,susceptible_color,susceptible_color,susceptible_color,susceptible_color, susceptible_color, susceptible_color
        , susceptible_color]
    nx.draw_networkx_nodes(G, pos=pos1, ax=axs[0,0], node_color=node_colors_step1)
    nx.draw_networkx_edges(G, pos=pos1, ax=axs[0,0], edge_color='black')
    nx.draw_networkx_labels(G, pos=pos1, ax=axs[0,0], labels=node_labels)
    axs[0,0].axis('off')

    node_colors_step1 = [infected_color, infected_color, infected_color, infected_color, susceptible_color,
                         susceptible_color, susceptible_color, susceptible_color, susceptible_color, susceptible_color,
                         susceptible_color, susceptible_color]
    node_labels = {0: '0', 1: '1', 2: '1', 3: '1', 4: '', 5: '', 6: '', 7: '', 8: '', 9:'', 10:'', 11:''}
    nx.draw_networkx_nodes(G, pos=pos1, ax=axs[0,1], node_color=node_colors_step1)
    nx.draw_networkx_edges(G, pos=pos1, ax=axs[0,1], edge_color='black')
    nx.draw_networkx_labels(G, pos=pos1, ax=axs[0,1], labels=node_labels)
    axs[0,1].axis('off')

    node_colors_step1 = [infected_color, infected_color, infected_color, infected_color, infected_color,
                         infected_color, susceptible_color, infected_color, susceptible_color, susceptible_color,
                         susceptible_color, susceptible_color]
    node_labels = {0: '0', 1: '1', 2: '1', 3: '1', 4: '2', 5: '2', 6: '', 7: '2', 8: '', 9:'', 10:'', 11:''}
    nx.draw_networkx_nodes(G, pos=pos1, ax=axs[1,0], node_color=node_colors_step1)
    nx.draw_networkx_edges(G, pos=pos1, ax=axs[1,0], edge_color='black')
    nx.draw_networkx_labels(G, pos=pos1, ax=axs[1,0], labels=node_labels)
    axs[1,0].axis('off')

    node_colors_step1 = [infected_color, infected_color, infected_color, infected_color, infected_color,
                         infected_color, susceptible_color, infected_color, infected_color, susceptible_color,
                         susceptible_color, susceptible_color]
    node_labels = {0: '0', 1: '1', 2: '1', 3: '1', 4: '2', 5: '2', 6: '', 7: '2', 8: '3', 9:'', 10:'', 11:''}
    nx.draw_networkx_nodes(G, pos=pos1, ax=axs[1,1], node_color=node_colors_step1)
    nx.draw_networkx_edges(G, pos=pos1, ax=axs[1,1], edge_color='black')
    nx.draw_networkx_labels(G, pos=pos1, ax=axs[1,1], labels=node_labels)

    axs[1,1].axis('off')

    plt.savefig('cartoon_network.svg', fmt='svg')
    plt.show()



