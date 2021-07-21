import matplotlib.pyplot as plt
from matplotlib import rc
from src import degree_distributions
from src import plotting_util
from src import figs_for_paper
from analysis import covid
import numpy as np
import random
from src import pgf_formalism

if __name__ == '__main__':
    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)
    print('pgfs yay!')
    covid.contour_fig3()
    plt.semilogy()
    plt.legend(loc='upper left')
    plt.show()

    #TODO do the fourier transform version with the intervention ones
    # power_law_q2 = degree_distributions.power_law_degree_distrb(2000, mu=10)
    # power_law_q2 = degree_distributions.binomial_degree_distb(500, 2.5)
    # pgf_formalism.compute_phase_space(21, 400, power_law_q2, 0.8, True, [2, 3, 4, 5, 6, 10, 12, 15, 20], '../data/talk_results/plaw_mod_random_vax{0}',
    #                          intervention_gen=4, vacc_pop=0.5, intervention_trans=0.0, do_non_interv=True, do_interv=True, intervention_type="random")
    # plotting_util.outbreak_size_curves([2, 4, 6, 10, 15], 200, '../data/talk_results/plaw_mod_random_vax{0}.txt', '../data/talk_results/plaw_mod_random_vax{0}_intv.txt', same_plot=True)
    # plt.show()
    #### TESTING THINGS

    # Old plots w data for NERCCS talk, can change colors etc if need be
    # plotting_util.plots_for_nerccs_talk([3, 6, 10, 15], 400, '../data/talk_results/pl_T8_no_intervene.txt',
    #                                     '../data/talk_results/powerlawrollout456{0}.txt',
    #                                     '../data/talk_results/universal_pl_T8_gen4_T4_intervene.txt',
    #                                     '../data/talk_results/powerlawuniversal4_T4{0}_intv.txt',
    #                                     same_plot=True, normalize_axis_x=False,
    #                                     inset_label=f'$T_g \\rightarrow T_4=0.4$')
    # plt.show()

    # Running simulations for interventions
    figs_for_paper.random_vax_exp('ensuring_software', num_sims=20000, num_nodes=10000)
    print("done")
    # ## Random Vaccination
    #TBD
    plotting_util.plots_for_nerccs_talk([3, 5, 6, 15], 400,  '../data/paper/plaw_T8_10k_120ksims_q3_generational.txt', '../data/talk_results/plaw_mod_random_vax{0}.txt',
                                                    '../data/testing/plaw_random_vax_10knodes_20k_sims_intervene.txt', '../data/talk_results/plaw_mod_random_vax{0}_intv.txt',
                                                    same_plot=True, normalize_axis_x=False, inset_label=f'$T_g \\rightarrow T_4=0$ for \n 50\% Vaccinated')
    plt.show()


    ###### done

    # covid.get_custom_cmap()
    # covid.contour_figs()
    # figs_for_paper.network_drawing()
    covid.contour_fig3()
    # covid.contour_fig2()

    # covid.exploratory_state_plot('Arkansas')
    # covid.exploratory_state_plot('Arizona') #use
    # covid.exploratory_state_plot('South Dakota') #use
    # covid.exploratory_state_plot('Delaware')
    # covid.exploratory_state_plot('New Jersey')
    # covid.exploratory_state_plot('Hawaii') #use
    # covid.exploratory_state_plot('Kentucky')
    # covid.exploratory_state_plot('Massachusetts') #use
    plt.semilogy()
    plt.legend(loc='upper left')
    plt.show()
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
                                                    plot_distribution_inset=False, legend_on=False)
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