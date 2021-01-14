import numpy as np
import matplotlib.pyplot as plt
import matplotlib as m
import matplotlib.colors as colors
import matplotlib.patches as mpatches


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
    # ax.imshow(psi_g[:60][:, :100], cmap=cmap, norm=colors.PowerNorm(gamma=0.05, vmin=0, vmax=max(psi_g[0])), label='$gen='+str(gen)+'$')  # gen 5
    ax.imshow(psi_g[:60][:, :100], cmap=cmap, norm=plt.cm.colors.SymLogNorm(linthresh=0.00005, vmax=0.4, vmin=0.000),
              label='$gen=' + str(gen) + '$')  # gen 5
    red_patch = mpatches.Patch(color='white', alpha=0.001, label='$gen=' + str(gen) + '$')
    plt.legend(handles=[red_patch], loc='upper right')
    # ax.imshow(psi_g[:80][:, :150], cmap=cm)  # gen 5
    ax.invert_yaxis()
    # plt.title('Phase Space at Generation '+str(gen)+' of '+str(title_label))
    plt.ylabel('$m$', fontsize=16)
    plt.xlabel('$s$', fontsize=16)
    # plt.legend(loc='upper right')
    # plt.savefig(str(title_label)+str(gen)+'_phase_space.png')
    plt.show()


def phaseSpace_from_data(fname, gen, plot_title):
    psi_g = np.loadtxt(fname, delimiter=',')
    inverted_s_m = psi_g.T
    plot_psi(inverted_s_m, gen, plot_title)


def plot_sims_vs_analytical_multigens(list_of_gens, x_lim, fname_sim_results, fname_predict_format, fname_sim_results_int=None, fname_predict_format_int=None, same_plot=False):
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
        fname_sims = fname_sim_results
        if fname_sim_results_int is not None and fname_predict_format_int is not None:
            fname_predict_interv = fname_predict_format_int.format(gen)
            fname_sims_interv = fname_sim_results_int
        else:
            fname_sims_interv = None
            fname_predict_interv = None
        plot_sims_vs_analytical_outbreak_sizes(gen, x_lim, fname_sims, fname_predict, fname_sims_interv,
                                               fname_predict_interv, color_key, same_plot)
    if same_plot:
        plt.show()


def plot_sims_vs_analytical_outbreak_sizes(gen, x_lim, fname_sim_results, fname_predict, fname_sim_results_interv,
                                           fname_predict_interv, color_key, same_plot=False):
    # Rough method for plotting simulations vs analytical probabilities of outbreak size.
    # Modify as needed for existing files or re-generation of probability results
    plot_intervention = False
    if fname_sim_results_interv is not None and fname_predict_interv is not None:
        plot_intervention = True

    data = np.loadtxt(fname_sim_results, delimiter=',')
    time_series = data[gen][2:x_lim]
    plt.plot(np.arange(2, x_lim), time_series, color=color_key[gen], alpha=0.95, ls='-', lw=.6)

    if plot_intervention:
        data_int = np.loadtxt(fname_sim_results_interv, delimiter=',')
        time_series_int = data_int[gen][2:x_lim]
        plt.plot(np.arange(2, x_lim), time_series_int, color=color_key[gen], ls='--', alpha=0.95, lw=.6)

    plt.semilogy()
    plt.ylim(.0001, .1)
    # plt.show()

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
