import matplotlib as m
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

import src.pgf_formalism
from src import degree_distributions


### Below are three sets of plotting functions; tools for plotting simulation results, tools for plotting
### analytical results (phase space and cumulative outbreak sizes per generation) and tools for comparing analytical and
### simulated results.

### SIMULATION PLOTTING TOOLS
def graph_infection_size_distribution_by_gen(list_of_gens, x_lim, filename,
                                             gen_emergence_times=None,
                                             intervention_filename=None):
    intervention_comparison_true = False
    if intervention_filename is not None:
        intervention_comparison_true = True
    color_key = {}
    colors = ['red', 'orange', 'green', 'blue', 'goldenrod', 'purple', 'pink', 'teal', 'black', 'gold', 'chocolate',
              'dodgerblue', 'darkslategray', 'mediumorchid', 'magenta', 'tomato', 'midnightblue',
              'cadetblue', 'crimson']
    colors = ['black', 'blue', 'red', 'purple', 'orange', 'brown']
    for i in range(len(list_of_gens)):
        gen = list_of_gens[i]
        color = colors[0]
        colors.remove(color)
        color_key[gen] = color

    data_no_intervention = np.loadtxt(filename, delimiter=',')

    for gen in list_of_gens:
        time_label=''
        if gen_emergence_times is not None:
            time_label = int(gen_emergence_times[gen])
        time_series = data_no_intervention[gen][2:x_lim]
        plt.plot(time_series, label=f'Cumulative infections up to birth of gen$=${str(gen+1)}', color=color_key[gen], alpha=0.7, lw=1)

    if intervention_comparison_true:
        data_intervention = np.loadtxt(intervention_filename, delimiter=',')
        for gen in list_of_gens:
            time_series_int = data_intervention[gen][2:x_lim]
            plt.plot(time_series_int, color=color_key[gen], alpha=0.75, ls='--', lw=1)

    plt.legend(loc='upper right')
    plt.xlabel('number infected $s$ at generation $g$')
    plt.ylabel('$p_s^g$')
    # plt.title(
    #     'Effects of simulations with intervention from $T=.2$ to $T=.1$ at $g=3$ on Binomial Degree Distribution Network')
    plt.semilogy()
    # plt.ylim(.0001, .1)
    # plt.title('Created from saved data')
    # plt.savefig('p_s_g_distribution_intervention.png')
    # plt.show()

    # TODO add an inset plot option with total number of infections?


#### PHASE SPACE PLOTTING TOOLS
def plot_psi(psi_g, gen, title_label):
    cmap = plt.cm.hot(np.linspace(1, 0, 100000))
    cmap = m.colors.ListedColormap(cmap[:, :-1])

    fig, ax = plt.subplots()
    ax.imshow(psi_g[:50][:, :80], cmap=cmap, norm=plt.cm.colors.SymLogNorm(linthresh=0.00005, vmax=0.4, vmin=0.000),
              label='$gen=' + str(gen) + '$')  # gen 5
    red_patch = mpatches.Patch(color='white', alpha=0.001, label='$gen=' + str(gen) + '$')
    # plt.legend(handles=[red_patch], loc='upper right')
    ax.invert_yaxis()
    # plt.title('Phase Space at Generation '+str(gen)+' of '+str(title_label))
    plt.ylabel('$m$ new ', fontsize=20)
    plt.xlabel('$s$ cumulative', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.tight_layout()
    plt.show()

### For NERCCS talk
def colorlist2(c1, c2, num):
    l = np.linspace(0, 1, num)
    a = np.abs(np.array(c1) - np.array(c2))
    m = np.min([c1, c2], axis=0)
    s = np.sign(np.array(c2) - np.array(c1)).astype(int)
    s[s == 0] = 1
    r = np.sqrt(np.c_[(l * a[0] + m[0])[::s[0]], (l * a[1] + m[1])[::s[1]], (l * a[2] + m[2])[::s[2]]])
    return r

def plot_psi_compressed(psi_g, gen):
    cmap = plt.cm.hot(np.linspace(1, 0, 100000))[::-1]
    cmap = m.colors.ListedColormap(cmap[:, :])

    compressed = np.sum(psi_g, axis=0)[2:60]
    x = np.linspace(0, 2 * np.pi, len(compressed))
    x = np.arange(2,60)
    # x_labels = np.arange(0, 60, 60/(len(x)))
    # plt.xticks(x, x_labels)
    y = compressed

    # cmap = LinearSegmentedColormap.from_list("", colorlist2((1, 0, 0), (0, 0, 1), 100))

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-2], points[1:-1], points[2:]], axis=1)

    lc = LineCollection(segments, cmap=cmap, linewidth=10)
    lc.set_array(x)
    plt.gca().add_collection(lc)
    plt.gca().autoscale()
    plt.ylim(.001, .1)
    plt.semilogy()
    plt.xlabel('$s$ cumulative', fontsize=20)
    plt.ylabel('Probability', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.tight_layout()
    plt.show()

    fig, ax = plt.subplots()
    # plt.plot(np.sum(psi_g, axis=0)[2:50], cmap=cmap, norm=plt.cm.colors.SymLogNorm(linthresh=0.00005, vmax=0.4, vmin=0.000))
    # plt.semilogy()
    # plt.show()

def plot_extinct_prob(psi, g_list = (3,5,7,10), sum_ax = 2, x_ax="x axis", y_ax="y axis", hide_zero = False):
    s_list = range(1,psi.shape[1])
    Psi = psi.sum(axis=sum_ax)
    if hide_zero:
        Psi = np.where(Psi==0, -1, Psi)
    plt.figure()
    markers = ['ko', 'bs', 'rD', 'mv']
    for g in g_list:
        lbl = "gen" + str(g)
        plt.plot(s_list, np.log(Psi[g][s_list[0]:s_list[-1]+1]), markers[g_list.index(g)], label=lbl)
    plt.xlabel(x_ax, fontsize=20)
    plt.ylabel(y_ax, fontsize=20)
    plt.xticks(range(0,len(s_list),50),fontsize=18)
    #plt.ylim((0,np.max(Psi)+0.05))
    plt.yticks(fontsize=18)
    plt.legend()
    plt.tight_layout()
    plt.show() 

def phaseSpace_from_data(fname, gen, plot_title):
    psi_g = np.loadtxt(fname, delimiter=',')
    inverted_s_m = psi_g.T
    plot_psi(inverted_s_m, gen, plot_title)
    plot_psi_compressed(inverted_s_m, gen)


def distribution_heatmap(num_gens, s_lim, degree_distribution, transmissibility):
    # This is both a computational and visualization function
    # Since the full matrix results for all generations from 0 to num_gens is required, and cumbersome to save a single
    # file for each generation, the entire phase space is computed here in memory and plotted.
    heatmap_m = np.zeros((num_gens, s_lim))
    heatmap_s = np.zeros((num_gens, s_lim))

    initProb = 1
    # If desired, specify intervention parameters below to obtain the phase space with intervention:
    all_psi_results = src.pgf_formalism.Psi(degree_distribution, initProb, num_gens, s_lim, s_lim, transmissibility, 4,
                                            0.001)

    for gen in range(1, num_gens):
        inverted_s_m = all_psi_results[gen].T
        s_marginal = np.sum(inverted_s_m, axis=0)
        m_marginal = np.sum(inverted_s_m, axis=1)
        s_marginal = s_marginal / np.sum(s_marginal)  # normalize
        m_marginal = m_marginal / np.sum(m_marginal)  # normalize
        heatmap_m[gen] = m_marginal
        heatmap_s[gen] = s_marginal

    # FIRST FIGURE:
    # x-axis: generations, y-axis: m, number infected during that generation
    # Heat map encodes probability distributions over m for each generation g, as a verticle bar from yellow to red
    cmap = plt.cm.hot(np.linspace(1, 0, 100000))
    cmap = m.colors.ListedColormap(cmap[:, :-1])

    fig, ax = plt.subplots()
    ax.axis('equal')
    ax.imshow(heatmap_m.T, cmap=cmap, norm=plt.cm.colors.SymLogNorm(linthresh=0.00005, vmax=0.4, vmin=0.000),
              aspect='auto')
    ax.invert_yaxis()
    plt.ylabel('$m$', fontsize=16)
    plt.xlabel('$gen$', fontsize=16)
    plt.xticks(np.arange(0, num_gens, 10), np.arange(0, num_gens, 10))
    plt.show()

    # SECOND FIGURE:
    # x-axis: generations, y-axis: s, total number infected BY that generation
    # Heat map encodes probability distributions over s for each generation g, as a verticle bar from yellow to red
    fig, ax = plt.subplots()
    ax.axis('equal')
    ax.imshow(heatmap_s.T, cmap=cmap, norm=plt.cm.colors.SymLogNorm(linthresh=0.00005, vmax=0.4, vmin=0.000),
              aspect='auto')
    ax.invert_yaxis()
    plt.ylabel('$s$', fontsize=16)
    plt.xlabel('$gen$', fontsize=16)
    plt.xticks(np.arange(0, num_gens, 10), np.arange(0, num_gens, 10))
    plt.show()


def outbreak_size_curves(list_of_gens, x_lim, fname_predict_format,
                         fname_predict_format_interv, same_plot=False):
    print('Analytical probability of total number infectives s at generations g with and without intervention')
    # Method for s-slices of the total phase space
    color_key = {}
    colors = ['red', 'orange', 'green', 'blue', 'purple', 'teal', 'black', 'gold', 'chocolate',
              'dodgerblue', 'darkslategray', 'mediumorchid', 'magenta', 'tomato', 'midnightblue',
              'cadetblue', 'crimson']

    for i in range(len(list_of_gens)):
        gen = list_of_gens[i]
        color = np.random.choice(colors)
        colors.remove(color)
        color_key[gen] = color

    for gen in list_of_gens:
        fname_predict = fname_predict_format.format(gen)
        if fname_predict_format_interv is not None:
            fname_predict_interv = fname_predict_format_interv.format(gen)
        else:
            fname_predict_interv = None
        plot_intervention = False
        if fname_predict_interv is not None:
            plot_intervention = True

        psi_g = np.loadtxt(fname_predict, delimiter=',')
        inverted_s_m = psi_g.T
        ps_g_analytical = np.sum(inverted_s_m, axis=0)
        ps_g_analytical = ps_g_analytical / np.sum(ps_g_analytical)  # normalize
        label = '$g=' + str(gen) + '$'
        color = color_key[gen]
        plt.plot(np.arange(2, x_lim), ps_g_analytical[2:x_lim], label=label, color=color, linestyle='-', lw=.8)

        if plot_intervention:
            psi_g_int = np.loadtxt(fname_predict_interv, delimiter=',')
            inverted_s_m_int = psi_g_int.T
            ps_g_analytical_int = np.sum(inverted_s_m_int, axis=0)
            ps_g_analytical_int = ps_g_analytical_int / np.sum(ps_g_analytical_int)  # normalize
            label = '$g=' + str(gen) + '$'
            color = color_key[gen]
            plt.plot(np.arange(2, x_lim), ps_g_analytical_int[2:x_lim], label=label + ' intervention', color=color,
                     linestyle='--', lw=.8)

        plt.semilogy()
        plt.ylim(.0001, .1)
        plt.rcParams.update({'font.size': 12})
        plt.legend(loc='upper right')
        plt.xlabel('$s$- number nodes infected at generation $g$', fontsize=12)
        plt.ylabel('$p_s^g$', fontsize=12)
        # plt.rcParams.update({'font.size': 12})
        plt.title('Effects of Intervention', fontsize=10)
        if not same_plot:
            plt.show()

    if same_plot:
        plt.show()
    print('done')


# todo make this so you can add things to a plot, one method in order to plot analytical v sims or just sims

### SIMULATION WITH ANALYTICAL PLOTTING TOOLS
def plot_sims_vs_analytical_multigens(list_of_gens, x_lim, fname_sim_results, fname_predict_format,
                                      fname_sim_results_int=None, fname_predict_format_int=None, same_plot=False,
                                      normalize_axis_x=False, plot_distribution_inset=False, inset_to_plot=None,
                                      inset_title=None, grayscale=False):
    color_key = {}
    color_key_sims = {}
    colors = ['red', 'orange', 'green', 'blue', 'purple', 'teal', 'black', 'gold', 'chocolate',
              'dodgerblue', 'darkslategray', 'mediumorchid', 'magenta', 'tomato', 'midnightblue',
              'cadetblue', 'crimson']

    cmap_blues = plt.cm.get_cmap('Blues')
    cmap_reds = plt.cm.get_cmap('Oranges')
    cmap_reds = plt.cm.get_cmap('brg')
    cmap_reds = plt.cm.get_cmap('CMRmap')
    cmap_reds = plt.cm.get_cmap('autumn')
    cmap_reds = plt.cm.get_cmap('gist_heat')
    # cmap_reds = plt.cm.get_cmap('gnuplot')
    hex_colorbrewer = ['#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e']
    hex_colorbrewer = ['#a1dab4','#41b6c4','#2c7fb8','#253494']
    standrd_colors = ['black', 'blue', 'red', 'purple']
    # standrd_colors = ["#48acf0","#594236","#6f584b", "#93a3bc"] # too light "#ccdde2",
    standrd_colors = ["#001524", "#15616d", "#ff7d00", "#78290f"] #  "#ffecd1", too light
    photocopy_friendly = ['#e66101','#fdb863','#b2abd2','#5e3c99']
    sequential_bluegreen = ['#a1dab4','#41b6c4','#2c7fb8','#253494'] # '#ffffcc', too light
    sequential_orange = ['#fecc5c','#fd8d3c','#f03b20','#bd0026'] # '#ffffb2', too light
    cmap = plt.get_cmap('bone')
    indices = np.linspace(0, cmap.N, len(list_of_gens)+1)
    my_colors = [cmap(int(i)) for i in indices]
    standrd_colors = my_colors[:-1][::-1] #Shifted to avoid yellow


    style_key = {}
    styles = ['-',':', '-.', '--', '-']

    for i in range(len(list_of_gens)):
        gen = list_of_gens[i]
        color = colors[0]
        colors.remove(color)
        color_key[gen] = color
        color_key[gen] = cmap_reds(1 - (0.3 + (i/(1.5*len(list_of_gens)))))
        color_key[gen] = standrd_colors[i]
        color_key_sims[gen] = cmap_blues(0.3 + (i/(1.5*len(list_of_gens))))
        style = styles[0]
        styles.remove(style)
        style_key[gen] = style

    # fig, ax1 = plt.subplots(figsize=(14, 7))
    fig, ax1 = plt.subplots(figsize=(8, 6))

    for gen in list_of_gens:
        fname_predict = fname_predict_format.format(gen)
        fname_sims = fname_sim_results
        if fname_sim_results_int is not None and fname_predict_format_int is not None:
            fname_predict_interv = fname_predict_format_int.format(gen)
            fname_sims_interv = fname_sim_results_int
        else:
            fname_sims_interv = None
            fname_predict_interv = None
        plot_sims_vs_analytical_outbreak_sizes(fig, ax1, gen, x_lim, fname_sims, fname_predict, fname_sims_interv,
                                               fname_predict_interv, color_key, style_key, same_plot, normalize_axis_x,
                                               # grayscale, color_key_sims)
                                               grayscale)

    if plot_distribution_inset:
        right, bottom, width, height = [0.4, 0.65, 0.25, 0.3]
        ax2 = fig.add_axes([right, bottom, width, height])
        # FOR POWER LAW with mu=10:
        # power_law_dd = degree_distributions.power_law_degree_distrb(400, mu=10) #q=3
        ax2.plot(inset_to_plot[:15], color='black', label=inset_title)
        # FOR ERDOS-RENYI with k=2.5:
        # erdos_renyi = degree_distributions.binomial_degree_distb(400, 2.5) #q=k=2.5
        # ax2.plot(erdos_renyi[:15], color='black', label='$p_k \\approx \\frac{\\lambda^ke^{-\\lambda}}{k!}, \\lambda=2.5$')
        ax2.set_xlim(0, 14)
        ax2.set_xlabel('Degree $k$', fontsize=20)
        ax2.set_ylabel('Fraction of nodes', fontsize=20)
        # ax2.set_yticks(np.arange(0, 1, 0.25))
        ax2.set_xticks([0, 1, 2, 3, 5, 10])
        ax2.semilogy()
        ax2.legend(loc='upper right', fontsize=16, frameon=False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['bottom'].set_visible(True)
        ax2.spines['left'].set_visible(True)

    # TODO adjust for this particular figure
    # ax1.text(0.0005, 0.001, '$g=2$')
    # ax1.text(0.0019, 0.0008, '$g=5$')
    # ax1.text(0.008, 0.0006, '$g=10$')
    # ax1.text(0.009, 0.0005, '$g=20$')
    intervention_color = '0.5'
    intervention_color = 'cornflowerblue'
    regular_color = '0.1'
    regular_color = 'navy'
    regular_patch = mpatches.Patch(color=regular_color, label='No intervention')
    intervention_patch = mpatches.Patch(color=intervention_color, label='Intervention applied at gen 4')
    # extra_legend = ax1.legend(handles=[regular_patch, intervention_patch], loc=(.15, .85), frameon=False)

    # TODO these legend elements are very particular to the static figure we made, change as needed
    if grayscale:
        legend_elements = [Line2D([0], [0], color='navy', lw=1, ls=style_key[2], label='Infections up to gen 2'),
                           Line2D([0], [0], color='navy', lw=1, ls=style_key[4], label='Infections up to gen 4'),
                           Line2D([0], [0], color='navy', lw=1, ls=style_key[10], label='Infections up to gen 10'),
                           regular_patch, intervention_patch]
        ax1.legend(handles=legend_elements, loc=(.1, .59), frameon=False, fontsize=18)
    plt.tight_layout()
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    ax1.tick_params(axis='y', labelrotation=0, labelsize=16)
    ax1.tick_params(axis='x', labelrotation=0, labelsize=16)
    ax1.legend(loc='upper right', fontsize=16, frameon=False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(True)
    ax1.spines['left'].set_visible(True)
    if same_plot:
        plt.legend(loc='upper right')
        plt.show()


def plot_sims_vs_analytical_outbreak_sizes(fig, ax1, gen, x_lim, fname_sim_results, fname_predict,
                                           fname_sim_results_interv,
                                           fname_predict_interv, color_key, style_key, same_plot=False,
                                           normalize_axis_x=False, grayscale=False, color_key_sims=None):
    # Rough method for plotting simulations vs analytical probabilities of outbreak size.
    # Modify as needed for existing files or re-generation of probability results

    intervention_color = '0.5'
    intervention_color = 'cornflowerblue'
    regular_color = '0.1'
    regular_color = 'navy'

    ax1.set_xlim(0, x_lim)

    x_start = 2

    plot_intervention = False
    if fname_sim_results_interv is not None and fname_predict_interv is not None:
        plot_intervention = True

    data = np.loadtxt(fname_sim_results, delimiter=',')
    time_series = data[gen][x_start:x_lim]
    x_vals = np.arange(x_start, x_lim)
    x_ticks = np.arange(2, x_lim, int((x_lim - 2) / 10))
    if normalize_axis_x:
        x_vals = x_vals / 10000
        x_ticks = x_ticks / 10000
        # ax1.set_xlim(2/10000, max(x_vals))
        ax1.set_xlim(0, max(x_vals))
    x_tick_labels = []
    for x in x_ticks:
        percent = np.round(x * 100, 2)
        str_percent = '${0}\\%$'.format(percent)
        x_tick_labels.append(str_percent)
    color = color_key[gen]
    if color_key_sims is not None:
        color = color_key_sims[gen]
    if grayscale:
        color = regular_color
    ax1.plot(x_vals, time_series, color=color, ls=style_key[gen], lw=1) #, label=f'gen {gen}')

    if plot_intervention:
        data_int = np.loadtxt(fname_sim_results_interv, delimiter=',')
        time_series_int = data_int[gen][x_start:x_lim]
        x_vals = np.arange(x_start, x_lim)
        if normalize_axis_x:
            x_vals = x_vals / 10000
        color = color_key[gen]
        alpha = 0.6
        if grayscale:
            color = intervention_color
            alpha = 1
        ax1.plot(x_vals, time_series_int, color=color, ls=style_key[gen], lw=1, alpha=alpha)

    ax1.semilogy()
    ax1.set_ylim(.00005, .1)
    if normalize_axis_x:
        plt.xticks(x_ticks, x_tick_labels)
    # if normalize_axis_x:
    #     plt.xlim(0, 1)
    # plt.show()

    psi_g = np.loadtxt(fname_predict, delimiter=',')
    inverted_s_m = psi_g.T
    ps_g_analytical = np.sum(inverted_s_m, axis=0)
    ps_g_analytical = ps_g_analytical / np.sum(ps_g_analytical)  # normalize
    x_vals = np.arange(x_start, x_lim)
    if normalize_axis_x:
        x_vals = x_vals / 10000
        # TODO fix labeling issue
    label = 'Infections up to gen {0}'.format(gen)
    label = 's at generation {0}'.format(gen)
    color = color_key[gen]
    if grayscale:
        color = regular_color
    ax1.plot(x_vals, ps_g_analytical[x_start:x_lim], label=label, color=color, lw=2, ls=style_key[gen])

    if plot_intervention:
        psi_g_int = np.loadtxt(fname_predict_interv, delimiter=',')
        inverted_s_m_int = psi_g_int.T
        ps_g_analytical_int = np.sum(inverted_s_m_int, axis=0)
        ps_g_analytical_int = ps_g_analytical_int / np.sum(ps_g_analytical_int)  # normalize
        x_vals = np.arange(x_start, x_lim)
        if normalize_axis_x:
            x_vals = x_vals / 10000
        # label = '$g=' + str(gen) + '$'
        color = color_key[gen]
        alpha = 0.6
        if grayscale:
            color = intervention_color
            alpha = 1
        ax1.plot(x_vals, ps_g_analytical_int[x_start:x_lim], color=color,
                 ls=style_key[gen], lw=1, alpha=alpha)

    plt.rcParams.update({'font.size': 12})
    # ax1.add_artist(plt.legend(loc=(.002, .005)))
    # ax1.add_artist(plt.legend(loc=(.15, .69), frameon=False, facecolor='white', framealpha=1))
    # plt.legend(loc=(.15, .69), frameon=False, facecolor='white', framealpha=1)
    # TODO remove boxes around legends
    # regular_patch = mpatches.Patch(color=regular_color,  label='No intervention')
    # intervention_patch = mpatches.Patch(color=intervention_color,  label='Intervention applied at gen 4')
    # # extra_legend = ax1.legend(handles=[regular_patch, intervention_patch], loc=(.15, .85), frameon=False)
    #
    # legend_elements = [Line2D([0], [0], color='black', lw=1, ls=style_key[gen], label=label),
    #                    regular_patch, intervention_patch]
    # ax1.legend(handles=legend_elements, loc=(.15, .69))
    # extra_legend = ax1.legend(handles=[regular_patch, intervention_patch], loc=(.1, .85))
    # ax1.add_artist(extra_legend)
    ax1.set_xlabel('Cumulative infections $s$', fontsize=20)
    ax1.set_ylabel('Probability', fontsize=20)
    # plt.rcParams.update({'font.size': 12})
    # plt.title('Effects of Intervention', fontsize=10)
    print(f'gen {gen}')
    print(np.sum((ps_g_analytical[:400]))-(np.sum(data[gen][:400])))
    if not same_plot:
        plt.show()

def gen_vs_time_correlation(gen_distributions, time_distributions):
    max_len = min(len(gen_distributions), len(time_distributions))
    corr_coef_series = np.zeros(max_len)
    for gen in range(max_len):
        gen_dist = gen_distributions[gen]
        time_dist = time_distributions[gen]
        # x_lim = 50
        # plt.plot(gen_dist[:max_len])
        # plt.plot(time_dist[:max_len])
        # plt.semilogy()
        # plt.show()
        correlation = np.corrcoef(gen_dist, time_dist)
        corr_coef_series[gen] = correlation[0][1]
    return corr_coef_series



def plots_for_nerccs_talk(list_of_gens, x_lim, fname_sim_results, fname_predict_format,
                                      fname_sim_results_int=None, fname_predict_format_int=None, same_plot=False,
                                      normalize_axis_x=False, plot_distribution_inset=False, grayscale=False, inset_label=None):
    ## Plain axes:
    fig, ax1 = plt.subplots(figsize=(15, 7))
    ax1.set_xlim(0, x_lim)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    x_vals = np.arange(1, x_lim)
    x_ticks = np.arange(2, x_lim, int((x_lim - 2) / 10))
    ax1.set_xlim(0, max(x_vals))
    ax1.set_xticks(x_ticks)
    ax1.semilogy()
    ax1.set_ylim(.00005, .1)
    ax1.set_xlabel('Cumulative Infections', fontsize=20)
    ax1.set_ylabel('Probability', fontsize=20)
    plt.tight_layout()
    plt.show()

    ## Axes with theory and inset
    fig, ax1 = plt.subplots(figsize=(15, 7))
    ax1.set_xlim(0, x_lim)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    x_vals = np.arange(1, x_lim)
    x_ticks = np.arange(2, x_lim, int((x_lim - 2) / 10))
    ax1.set_xlim(0, max(x_vals))
    ax1.set_xticks(x_ticks)
    ax1.semilogy()
    ax1.set_ylim(.00005, .1)
    ax1.set_xlabel('Cumulative Infections', fontsize=20)
    ax1.set_ylabel('Probability', fontsize=20)
    if plot_distribution_inset:
        right, bottom, width, height = [0.4, 0.5, 0.25, 0.3]
        ax2 = fig.add_axes([right, bottom, width, height])
        power_law_dd = degree_distributions.power_law_degree_distrb(400)
        ax2.plot(power_law_dd[:15], color='black')
        ax2.set_xlim(0, 14)
        ax2.set_title('Network Degree Distribution', fontsize=20)
        ax2.set_xlabel('Degree $k$', fontsize=20)
        ax2.set_ylabel('Fraction of nodes', fontsize=20)
        # ax2.set_yticks(np.arange(0, 1, 0.25))
        ax2.set_xticks([0, 1, 2, 3, 5, 10])
        ax2.semilogy()
    plt.tight_layout()
    plt.show()

    ## Axes with just theory
    fig, ax1 = plt.subplots(figsize=(15, 7))
    ax1.set_xlim(0, x_lim)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    x_vals = np.arange(1, x_lim)
    x_ticks = np.arange(2, x_lim, int((x_lim - 2) / 10))
    ax1.set_xlim(0, max(x_vals))
    ax1.set_xticks(x_ticks)
    ax1.semilogy()
    ax1.set_ylim(.00005, .1)
    ax1.set_xlabel('Cumulative Infections', fontsize=20)
    ax1.set_ylabel('Probability', fontsize=20)

    color_key = {}
    colors = ['blue', 'teal', 'orange', 'crimson', 'purple', 'teal', 'black', 'gold', 'chocolate',
              'dodgerblue', 'darkslategray', 'mediumorchid', 'magenta', 'tomato', 'midnightblue',
              'cadetblue', 'crimson']
    colorcycle = sns.color_palette("mako_r", len(list_of_gens))

    ## Color test: Uncomment to test plot
    # colorcycle = sns.color_palette("mako_r", 4)
    # plt.plot(np.arange(100), np.sin(np.arange(100)) + 0, color=colorcycle[0])
    # plt.plot(np.arange(100), np.sin(np.arange(100)) + 1, color=colorcycle[1])
    # plt.plot(np.arange(100), np.sin(np.arange(100)) + 2, color=colorcycle[2])
    # plt.plot(np.arange(100), np.sin(np.arange(100)) + 3, color=colorcycle[3])
    # plt.tight_layout()
    # plt.show()

    style_key = {}
    styles = [':', '-.', '--', '-']

    for i in range(len(list_of_gens)):
        gen = list_of_gens[i]
        # color = colors[0]
        # colors.remove(color)
        color_key[gen] = colorcycle[i]
        style = styles[0]
        styles.remove(style)
        style_key[gen] = style

    for gen in list_of_gens:
        fname_predict = fname_predict_format.format(gen)
        fname_sims = fname_sim_results
        if fname_sim_results_int is not None and fname_predict_format_int is not None:
            fname_predict_interv = fname_predict_format_int.format(gen)
            fname_sims_interv = fname_sim_results_int
        else:
            fname_sims_interv = None
            fname_predict_interv = None
        nerccs_plot_sims_vs_analytical_outbreak_sizes(fig, ax1, gen, x_lim, None, fname_predict, None,
                                               None, color_key, style_key, same_plot, normalize_axis_x,
                                               grayscale)

    if plot_distribution_inset:
        right, bottom, width, height = [0.4, 0.5, 0.25, 0.3]
        ax2 = fig.add_axes([right, bottom, width, height])
        power_law_dd = degree_distributions.power_law_degree_distrb(400)
        ax2.plot(power_law_dd[:15], color='black')
        ax2.set_xlim(0, 14)
        ax2.set_title('Network Degree Distribution', fontsize=20)
        ax2.set_xlabel('Degree $k$', fontsize=20)
        ax2.set_ylabel('Fraction of nodes', fontsize=20)
        # ax2.set_yticks(np.arange(0, 1, 0.25))
        ax2.set_xticks([0, 1, 2, 3, 5, 10])
        ax2.semilogy()

    plt.tight_layout()
    plt.show()

    ## Axes with theory and sims
    fig, ax1 = plt.subplots(figsize=(15, 7))
    ax1.set_xlim(0, x_lim)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    x_vals = np.arange(1, x_lim)
    x_ticks = np.arange(2, x_lim, int((x_lim - 2) / 10))
    ax1.set_xlim(0, max(x_vals))
    ax1.set_xticks(x_ticks)
    ax1.semilogy()
    ax1.set_ylim(.00005, .1)
    ax1.set_xlabel('Cumulative Infections', fontsize=20)
    ax1.set_ylabel('Probability', fontsize=20)

    for gen in list_of_gens:
        fname_predict = fname_predict_format.format(gen)
        fname_sims = fname_sim_results
        if fname_sim_results_int is not None and fname_predict_format_int is not None:
            fname_predict_interv = fname_predict_format_int.format(gen)
            fname_sims_interv = fname_sim_results_int
        else:
            fname_sims_interv = None
            fname_predict_interv = None
        nerccs_plot_sims_vs_analytical_outbreak_sizes(fig, ax1, gen, x_lim, fname_sims, fname_predict, None,
                                               None, color_key, style_key, same_plot, normalize_axis_x,
                                               grayscale)
    plt.tight_layout()
    plt.show()

    ## Axes with theory and intervention theory
    fig, ax1 = plt.subplots(figsize=(15, 7))
    ax1.set_xlim(0, x_lim)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    x_vals = np.arange(1, x_lim)
    x_ticks = np.arange(2, x_lim, int((x_lim - 2) / 10))
    ax1.set_xlim(0, max(x_vals))
    ax1.set_xticks(x_ticks)
    ax1.semilogy()
    ax1.set_ylim(.00005, .1)
    ax1.set_xlabel('Cumulative Infections', fontsize=20)
    ax1.set_ylabel('Probability', fontsize=20)

    for gen in list_of_gens:
        fname_predict = fname_predict_format.format(gen)
        fname_sims = fname_sim_results
        if fname_sim_results_int is not None and fname_predict_format_int is not None:
            fname_predict_interv = fname_predict_format_int.format(gen)
            fname_sims_interv = fname_sim_results_int
        else:
            fname_sims_interv = None
            fname_predict_interv = None
        nerccs_plot_sims_vs_analytical_outbreak_sizes(fig, ax1, gen, x_lim, None, fname_predict, None,
                                               fname_predict_interv, color_key, style_key, same_plot, normalize_axis_x,
                                               grayscale)

    right, bottom, width, height = [0.4, 0.75, 0.25, 0.003]
    ax2 = fig.add_axes([right, bottom, width, height], frameon=False)
    ax2.set_xlabel(inset_label, fontsize=25)
    ax2.set_xticks([])
    ax2.set_yticks([])
    # power_law_dd = degree_distributions.power_law_degree_distrb(400)
    # ax2.plot(power_law_dd[:15], color='black')
    # ax2.set_xlim(0, 14)
    # ax2.set_xlabel('Degree $k$', fontsize=20)
    # ax2.set_ylabel('Fraction of nodes', fontsize=20)
    # ax2.set_yticks(np.arange(0, 1, 0.25))
    # ax2.set_xticks([0, 1, 2, 3, 5, 10])
    # ax2.semilogy()

    plt.tight_layout()
    plt.show()

    ## Axes with theory and intervention theory and intervention sims
    fig, ax1 = plt.subplots(figsize=(15, 7))
    ax1.set_xlim(0, x_lim)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    x_vals = np.arange(1, x_lim)
    x_ticks = np.arange(2, x_lim, int((x_lim - 2) / 10))
    ax1.set_xlim(0, max(x_vals))
    ax1.set_xticks(x_ticks)
    ax1.semilogy()
    ax1.set_ylim(.00005, .1)
    ax1.set_xlabel('Cumulative Infections', fontsize=20)
    ax1.set_ylabel('Probability', fontsize=20)

    for gen in list_of_gens:
        fname_predict = fname_predict_format.format(gen)
        fname_sims = fname_sim_results
        if fname_sim_results_int is not None and fname_predict_format_int is not None:
            fname_predict_interv = fname_predict_format_int.format(gen)
            fname_sims_interv = fname_sim_results_int
        else:
            fname_sims_interv = None
            fname_predict_interv = None
        nerccs_plot_sims_vs_analytical_outbreak_sizes(fig, ax1, gen, x_lim, None, fname_predict, fname_sims_interv,
                                               fname_predict_interv, color_key, style_key, same_plot, normalize_axis_x,
                                               grayscale)
    right, bottom, width, height = [0.4, 0.75, 0.25, 0.003]
    ax2 = fig.add_axes([right, bottom, width, height], frameon=False)
    ax2.set_xlabel(inset_label, fontsize=25)
    ax2.set_xticks([])
    ax2.set_yticks([])
    # power_law_dd = degree_distributions.power_law_degree_distrb(400)
    # ax2.plot(power_law_dd[:15], color='black')
    # ax2.set_xlim(0, 14)
    # ax2.set_xlabel('Degree $k$', fontsize=20)
    # ax2.set_ylabel('Fraction of nodes', fontsize=20)
    # ax2.set_yticks(np.arange(0, 1, 0.25))
    # ax2.set_xticks([0, 1, 2, 3, 5, 10])
    # ax2.semilogy()
    plt.tight_layout()
    plt.show()


def nerccs_plot_sims_vs_analytical_outbreak_sizes(fig, ax1, gen, x_lim, fname_sim_results, fname_predict,
                                           fname_sim_results_interv,
                                           fname_predict_interv, color_key, style_key, same_plot=False,
                                           normalize_axis_x=False, grayscale=False):
    # Rough method for plotting simulations vs analytical probabilities of outbreak size.
    # Modify as needed for existing files or re-generation of probability results

    intervention_color = '0.5'
    intervention_color = 'cornflowerblue'
    regular_color = '0.1'
    regular_color = 'navy'
    line_w=2

    ax1.set_xlim(0, x_lim)

    plot_intervention_predict = False
    plot_intervention_sims = False
    if fname_sim_results_interv is not None:
        plot_intervention_sims = True
    if fname_predict_interv is not None:
        plot_intervention_predict = True

    if fname_sim_results is not None:
        data = np.loadtxt(fname_sim_results, delimiter=',')
        time_series = data[gen][1:x_lim]
        x_vals = np.arange(1, x_lim)

        color = color_key[gen]

        ax1.plot(x_vals, time_series, color=color, ls='-', lw=line_w, alpha=0.4)

    if plot_intervention_sims:
        data_int = np.loadtxt(fname_sim_results_interv, delimiter=',')
        time_series_int = data_int[gen][1:x_lim]
        x_vals = np.arange(1, x_lim)

        color = color_key[gen]
        alpha = 0.4
        if grayscale:
            color = intervention_color
            alpha = 1
        ax1.plot(x_vals, time_series_int, color=color, ls='--', lw=line_w, alpha=alpha)

    psi_g = np.loadtxt(fname_predict, delimiter=',')
    inverted_s_m = psi_g.T
    ps_g_analytical = np.sum(inverted_s_m, axis=0)
    ps_g_analytical = ps_g_analytical / np.sum(ps_g_analytical)  # normalize
    x_vals = np.arange(1, x_lim)
    label = 'Infections by gen {0}'.format(gen)
    color = color_key[gen]
    if grayscale:
        color = regular_color
    ax1.plot(x_vals, ps_g_analytical[1:x_lim], label=label, color=color, ls='-', lw=line_w)

    if plot_intervention_predict:
        psi_g_int = np.loadtxt(fname_predict_interv, delimiter=',')
        inverted_s_m_int = psi_g_int.T
        ps_g_analytical_int = np.sum(inverted_s_m_int, axis=0)
        ps_g_analytical_int = ps_g_analytical_int / np.sum(ps_g_analytical_int)  # normalize
        x_vals = np.arange(1, x_lim)
        color = color_key[gen]
        alpha = 1.0
        ax1.plot(x_vals, ps_g_analytical_int[1:x_lim], color=color,
                 ls='--', lw=line_w, alpha=alpha)

    plt.rcParams.update({'font.size': 12})
    plt.legend(loc='upper right', frameon=False)
    if not same_plot:
        plt.tight_layout()
        plt.show()
