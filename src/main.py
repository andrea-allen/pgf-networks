import matplotlib.pyplot as plt
from matplotlib import rc
from src import degree_distributions
from src import plotting_util
from src import figs_for_paper
from analysis import covid
import numpy as np
from src import pgf_formalism

if __name__ == '__main__':
    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)
    print('pgfs yay!')

    # covid.contour_figs()
    # covid.contour_fig2()

    power_law_dd = degree_distributions.power_law_degree_distrb(400)
    poisson = degree_distributions.binomial_degree_distb(2000, 2.5)
    # figs_for_paper.combine_plaw_results()
    # figs_for_paper.power_law_exp('plaw_06_10_40knet', num_sims=200, num_nodes=40000, active_gen_sizes_on=True)
    figs_for_paper.results_plots('testing/plaw_06_10_40knet', q_degree=3.04, active_gen_sizes_on=True)
    plotting_util.plot_sims_vs_analytical_multigens([3, 4, 6, 10], 400, f'../data/paper/plaw_combo_125k_generational.txt',
                                                    '../data/paper/powerlaw_q3_T8_ifft_g{0}.txt',
                                                    inset_to_plot=power_law_dd, inset_title='$p_k = k^{-2}e^{-k/10}$',
                                                    same_plot=True, normalize_axis_x=False,
                                                    plot_distribution_inset=False)
    plt.show()

    plotting_util.plot_sims_vs_analytical_multigens([3, 4, 6, 10], 400, f'../data/testing/ER_05_25_generational.txt',
                                                    '../data/paper/ER_q2.5_T8_ifft_g{0}.txt', inset_to_plot=poisson,
                                                    inset_title='tbd$p_k = k^{-2}e^{-k/10}$',
                                                    same_plot=True, normalize_axis_x=False,
                                                    plot_distribution_inset=False)
    plt.show()






    # covid.contour_figs()
    covid.contour_fig2()

    power_law_dd = degree_distributions.power_law_degree_distrb(400)
    power_law_q3 = degree_distributions.power_law_degree_distrb(2000, mu=10) #q is 3, k is 1.7

    # figs_for_paper.power_law_exp('plaw_05_27_500sims', num_sims=500, num_nodes=10000, active_gen_sizes_on=True)
    figs_for_paper.results_plots(file_root='testing/plaw_05_27_500sims', q_degree=3.04, active_gen_sizes_on=True) #used just for the time plot because only need 10k results,not 120k

    poisson = degree_distributions.binomial_degree_distb(2000, 2.5)


    plotting_util.plot_sims_vs_analytical_multigens([3, 4, 6, 10], 400,  f'../data/testing/plaw_T8_10k_120ksims_q3_generational.txt',
                                                    '../data/testing/powerlaw_q3_T8_ifft_g{0}.txt', inset_to_plot=power_law_q3, inset_title='$p_k = k^{-2}e^{-k/10}$',
                                                    same_plot=True, normalize_axis_x=False, plot_distribution_inset=True)
    plt.show()
    # # #
    plotting_util.plot_sims_vs_analytical_multigens([15, 6, 8, 10], 400,  f'../data/testing/erdos_renyi_10k_60ksims_combo_generational.txt',
                                                    '../data/paper/ER_q2.5_T8_ifft_g{0}.txt', inset_to_plot=poisson, inset_title='$p_k \\approx \\frac{\\lambda^ke^{-\\lambda}}{k!}, \\lambda=2.5$',
                                                    same_plot=True, normalize_axis_x=False, plot_distribution_inset=True)
    plt.show()

    figs_for_paper.results_plots(file_root='paper/plaw_T8_10k_1ksims_q3_p7', q_degree=3.04, active_gen_sizes_on=True) #used just for the time plot because only need 10k results,not 120k