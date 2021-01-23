import numpy as np
import matplotlib.pyplot as plt
import matplotlib as m
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import src.pgf_formalism
from src import degree_distributions
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
from matplotlib.lines import Line2D

# TODO need to edit some of the captions and titles to be customizable/more descriptive

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
    # plt.title(
    #     'Effects of simulations with intervention from $T=.2$ to $T=.1$ at $g=3$ on Binomial Degree Distribution Network')
    plt.semilogy()
    plt.ylim(.0001, .1)
    # plt.title('Created from saved data')
    # plt.savefig('p_s_g_distribution_intervention.png')
    plt.show()

    # TODO add an inset plot option with total number of infections?


def plot_psi(psi_g, gen, title_label):
    cmap = plt.cm.hot(np.linspace(1, 0, 100000))
    cmap = m.colors.ListedColormap(cmap[:, :-1])

    fig, ax = plt.subplots()
    ax.imshow(psi_g[:60][:, :100], cmap=cmap, norm=plt.cm.colors.SymLogNorm(linthresh=0.00005, vmax=0.4, vmin=0.000),
              label='$gen=' + str(gen) + '$')  # gen 5
    red_patch = mpatches.Patch(color='white', alpha=0.001, label='$gen=' + str(gen) + '$')
    plt.legend(handles=[red_patch], loc='upper right')
    ax.invert_yaxis()
    plt.title('Phase Space at Generation '+str(gen)+' of '+str(title_label))
    plt.ylabel('$m$', fontsize=16)
    plt.xlabel('$s$', fontsize=16)
    plt.show()


def phaseSpace_from_data(fname, gen, plot_title):
    psi_g = np.loadtxt(fname, delimiter=',')
    inverted_s_m = psi_g.T
    plot_psi(inverted_s_m, gen, plot_title)

def distribution_heatmap(num_gens, s_lim, degree_distribution, transmissibility):
    # This is both a computational and visualization function
    # Since the full matrix results for all generations from 0 to num_gens is required, and cumbersome to save a single
    # file for each generation, the entire phase space is computed here in memory and plotted.
    heatmap_m = np.zeros((num_gens, s_lim))
    heatmap_s = np.zeros((num_gens, s_lim))

    initProb = 1
    # If desired, specify intervention parameters below to obtain the phase space with intervention:
    all_psi_results = src.pgf_formalism.Psi(degree_distribution, initProb, num_gens, s_lim, s_lim, transmissibility, 4, 0.001)

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


def plot_sims_vs_analytical_multigens(list_of_gens, x_lim, fname_sim_results, fname_predict_format,
                                      fname_sim_results_int=None, fname_predict_format_int=None, same_plot=False,
                                      normalize_axis_x=False, plot_distribution_inset=False, grayscale=False):
    color_key = {}
    colors = ['red', 'orange', 'green', 'blue', 'purple', 'teal', 'black', 'gold', 'chocolate',
              'dodgerblue', 'darkslategray', 'mediumorchid', 'magenta', 'tomato', 'midnightblue',
              'cadetblue', 'crimson']

    style_key = {}
    styles = [':', '-.', '--', '-']

    for i in range(len(list_of_gens)):
        gen = list_of_gens[i]
        color = np.random.choice(colors)
        colors.remove(color)
        color_key[gen] = color
        style = styles[0]
        styles.remove(style)
        style_key[gen] = style


    fig, ax1 = plt.subplots(figsize=(14,6))

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
                                               fname_predict_interv, color_key, style_key, same_plot, normalize_axis_x, grayscale)

    if plot_distribution_inset:
        right, bottom, width, height = [0.6, 0.4, 0.25, 0.3]
        ax2 = fig.add_axes([right, bottom, width, height])
        power_law_dd = degree_distributions.power_law_degree_distrb(10000)
        ax2.plot(power_law_dd[:15], color='black')
        ax2.set_xlim(0,14)
        ax2.set_xlabel('Degree $k$')
        ax2.set_ylabel('Fraction of nodes')
        # ax2.set_yticks(np.arange(0, 1, 0.25))
        ax2.set_xticks([0, 1, 2, 3, 5, 10])
        ax2.semilogy()

    # TODO adjust for this particular figure
    # ax1.text(0.0005, 0.001, '$g=2$')
    # ax1.text(0.0019, 0.0008, '$g=5$')
    # ax1.text(0.008, 0.0006, '$g=10$')
    # ax1.text(0.009, 0.0005, '$g=20$')
    intervention_color = '0.5'
    regular_color = '0.1'
    regular_patch = mpatches.Patch(color=regular_color,  label='No intervention')
    intervention_patch = mpatches.Patch(color=intervention_color,  label='Intervention applied at gen 4')
    # extra_legend = ax1.legend(handles=[regular_patch, intervention_patch], loc=(.15, .85), frameon=False)

    legend_elements = [Line2D([0], [0], color='black', lw=1, ls=style_key[2], label = 'Infections up to gen 2'),
                       Line2D([0], [0], color='black', lw=1, ls=style_key[4], label='Infections up to gen 4'),
                       Line2D([0], [0], color='black', lw=1, ls=style_key[10], label='Infections up to gen 10'),
                       regular_patch, intervention_patch]
    ax1.legend(handles=legend_elements, loc=(.1, .69), frameon=False)
    plt.tight_layout()
    if same_plot:
        plt.show()


def plot_sims_vs_analytical_outbreak_sizes(fig, ax1, gen, x_lim, fname_sim_results, fname_predict, fname_sim_results_interv,
                                           fname_predict_interv, color_key, style_key, same_plot=False, normalize_axis_x=False, grayscale=False):
    # Rough method for plotting simulations vs analytical probabilities of outbreak size.
    # Modify as needed for existing files or re-generation of probability results

    intervention_color = '0.5'
    regular_color = '0.1'

    ax1.set_xlim(0, x_lim)

    plot_intervention = False
    if fname_sim_results_interv is not None and fname_predict_interv is not None:
        plot_intervention = True

    data = np.loadtxt(fname_sim_results, delimiter=',')
    time_series = data[gen][1:x_lim]
    x_vals = np.arange(1, x_lim)
    x_ticks = np.arange(2, x_lim, int((x_lim-2)/10))
    if normalize_axis_x:
        x_vals = x_vals/10000
        x_ticks = x_ticks/10000
        # ax1.set_xlim(2/10000, max(x_vals))
        ax1.set_xlim(0, max(x_vals))
    x_tick_labels = []
    for x in x_ticks:
        percent = np.round(x*100, 2)
        str_percent = '${0}\\%$'.format(percent)
        x_tick_labels.append(str_percent)
    color = color_key[gen]
    if grayscale:
        color = regular_color
    ax1.plot(x_vals, time_series, color=color,ls=style_key[gen], lw=1)

    if plot_intervention:
        data_int = np.loadtxt(fname_sim_results_interv, delimiter=',')
        time_series_int = data_int[gen][1:x_lim]
        x_vals = np.arange(1, x_lim)
        if normalize_axis_x:
            x_vals = x_vals / 10000
        color = color_key[gen]
        if grayscale:
            color = intervention_color
        ax1.plot(x_vals, time_series_int, color=color, ls=style_key[gen], lw=1)

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
    x_vals = np.arange(1, x_lim)
    if normalize_axis_x:
        x_vals = x_vals / 10000
        #TODO fix labeling issue
    label = 'Infections up to gen {0}'.format(gen)
    color = color_key[gen]
    if grayscale:
        color = regular_color
    ax1.plot(x_vals, ps_g_analytical[1:x_lim], label=label, color=color, ls=style_key[gen], lw=1)

    if plot_intervention:
        psi_g_int = np.loadtxt(fname_predict_interv, delimiter=',')
        inverted_s_m_int = psi_g_int.T
        ps_g_analytical_int = np.sum(inverted_s_m_int, axis=0)
        ps_g_analytical_int = ps_g_analytical_int / np.sum(ps_g_analytical_int)  # normalize
        x_vals = np.arange(1, x_lim)
        if normalize_axis_x:
            x_vals = x_vals / 10000
        # label = '$g=' + str(gen) + '$'
        color = color_key[gen]
        if grayscale:
            color = intervention_color
        ax1.plot(x_vals, ps_g_analytical_int[1:x_lim], color=color,
                 ls=style_key[gen], lw=1)

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
    ax1.set_xlabel('Proportion of population cumulatively infected', fontsize=14)
    ax1.set_ylabel('Probability', fontsize=14)
    # plt.rcParams.update({'font.size': 12})
    # plt.title('Effects of Intervention', fontsize=10)
    if not same_plot:
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
