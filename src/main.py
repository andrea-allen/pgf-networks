import time

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy import stats

from analysis import ensemble
from src import degree_distributions
from src import pgf_formalism
from src import plotting_util
from src import figs_for_paper

import epintervene

if __name__ == '__main__':
    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)
    print('pgfs yay!')

    # figs_for_paper.combine_data()

    power_law_dd = degree_distributions.power_law_degree_distrb(400)
    power_law_q3 = degree_distributions.power_law_degree_distrb(2000, mu=10) #q is 3, k is 1.7
    poisson = degree_distributions.binomial_degree_distb(2000, 2.5)
    #
    # pgf_formalism.compute_phase_space(num_gens=17, num_nodes=2000, degree_distribution=power_law_q3, transmissibility=0.8,
    #                                   save_results=True, gens_to_save=[2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 16],
    #                                   file_fmt_to_save='../data/testing/powerlaw_q3_T8_ifft_g{0}',
    #                                     do_non_interv=True, do_interv=False)

    # plotting_util.plot_sims_vs_analytical_multigens([3, 6, 10, 12], 400,  f'../data/testing/plaw_T8_10k_80ksims_generational.txt',
    #                                                 '../data/testing/powerlaw_T8_g{0}.txt',
    #                                                 same_plot=True, normalize_axis_x=False)
    # plt.show()

    # plotting_util.plot_sims_vs_analytical_multigens([3, 4, 6, 10], 400,  f'../data/testing/plaw_T8_10k_120ksims_q3_generational.txt',
    #                                                 '../data/testing/powerlaw_q3_T8_ifft_g{0}.txt',
    #                                                 same_plot=True, normalize_axis_x=False, plot_distribution_inset=True)
    # plt.show()
    # #
    # plotting_util.plot_sims_vs_analytical_multigens([15, 6, 8, 10], 400,  f'../data/testing/erdos_renyi_10k_60ksims_combo_generational.txt',
    #                                                 '../data/testing/poiss_T8_appx800_g{0}.txt',
    #                                                 same_plot=True, normalize_axis_x=False, plot_distribution_inset=True)
    # plt.show()
    # plt.figure('correctedps')
    # plotting_util.plot_sims_vs_analytical_multigens([3, 6, 8, 10], 400,  f'../data/testing/erdos_renyi_10k_60ksims_combo_generational.txt',
    #                                                 '../data/testing/poiss_T8_appx2k_uptom_g{0}.txt',
    #                                                 same_plot=True, normalize_axis_x=False, plot_distribution_inset=True)
    # plt.show()

    # ensemble.run_ensemble_intervention_effects(power_law_dd, '../data/testing/pl_20k', 30000, 20000,
    #                                            init_T=0.8, intervention_type="none", run_regular=True)

    start_time = time.time()
    # figs_for_paper.erdos_renyi_exp(file_root='erdos_renyi_10k_10ksims_pt1', num_sims=10000, num_nodes=10000, kill_by=16)
    # figs_for_paper.erdos_renyi_exp(file_root='erdos_renyi_10k_10ksims_pt2', num_sims=10000, num_nodes=10000, kill_by=16)
    # figs_for_paper.erdos_renyi_exp(file_root='erdos_renyi_10k_10ksims_pt3', num_sims=10000, num_nodes=10000, kill_by=16)
    # figs_for_paper.erdos_renyi_exp(file_root='erdos_renyi_10k_10ksims_pt4', num_sims=10000, num_nodes=10000, kill_by=16)
    # figs_for_paper.erdos_renyi_exp(file_root='erdos_renyi_10k_10ksims_pt5', num_sims=10000, num_nodes=10000, kill_by=16)
    # figs_for_paper.erdos_renyi_exp(file_root='erdos_renyi_10k_10ksims_pt6', num_sims=10000, num_nodes=10000, kill_by=16)
    # figs_for_paper.erdos_renyi_exp(file_root='theory_testing_10k_2k_p2', num_sims=2000, num_nodes=10000, kill_by=16)
    # figs_for_paper.erdos_renyi_exp(file_root='theory_testing_10k_2k_p3', num_sims=2000, num_nodes=10000, kill_by=16)
    # figs_for_paper.power_law_exp(file_root='plaw_T8_10k_1ksims_q3_p7', num_sims=1000, num_nodes=10000, active_gen_sizes_on=True)
    # figs_for_paper.power_law_exp(file_root='plaw_T8_10k_40ksims_q3_p5', num_sims=40000, num_nodes=10000)
    # figs_for_paper.power_law_exp(file_root='plaw_T8_10k_40ksims_q3_p6', num_sims=40000, num_nodes=10000)
    # figs_for_paper.power_law_exp(file_root='plaw_T8_20k_40ksims_p2', num_sims=4000, num_nodes=20000)
    # figs_for_paper.tree_exp(num_nodes=8191, num_sims=300)
    print(f'took {time.time()-start_time} seconds')
    # figs_for_paper.results_plots(file_root='plaw_T8_10k_30ksims', q_degree=1.65)
    # figs_for_paper.results_plots(file_root='poiss_T8_5k_q_6_gamma1_g_over_b', q_degree=6)
    # figs_for_paper.results_plots(file_root='erdos_renyi_10k_10ksims_pt6', q_degree=2.5)
    # figs_for_paper.results_plots(file_root='plaw_T8_10k_40ksims_q3_p6', q_degree=3.04, active_gen_sizes_on=False)
    # figs_for_paper.results_plots(file_root='plaw_T8_10k_1ksims_q3_p7', q_degree=3.04, active_gen_sizes_on=True)
    figs_for_paper.results_plots(file_root='plaw_T8_10k_40ksims_q3_p6', q_degree=3.04, active_gen_sizes_on=False)
    # figs_for_paper.results_plots(file_root='plaw_T8_10k_30ksims_p2', q_degree=1.65)

    #### NICHOLAS
    # power_law_dd = degree_distributions.power_law_degree_distrb(400)
    # extnct_cake = pgf_formalism.compute_extinct_prob_all(power_law_dd, T=0.8)
    # plotting_util.plot_psi(extnct_cake[3], 3, 'gen 3 extct cake')
    # plt.show()

    ### END NICHOLAS


    # ensemble.plotting_exercise()
    print('done trying something')

    power_law_dd = degree_distributions.power_law_degree_distrb(400)
    power_law_dd10 = degree_distributions.power_law_degree_distrb(400, alpha=2, mu=30)
    power_law_dd100 = degree_distributions.power_law_degree_distrb(400, alpha=1)
    poisson = degree_distributions.binomial_degree_distb(400, 2.5)
    crazy_line_dd = degree_distributions.chain_degree_dist(10)

    #### config model:
    degree_seq = np.array(
        [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 1, 1, 1, 1, 1, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1,
         1, 1])
    np.mean([degree_seq[k] - 1 for k in range(len(degree_seq))])
    distrb = np.zeros(10)
    for i in range(len(degree_seq)):
        distrb[degree_seq[i]] +=1

    distrb = distrb/len(degree_seq)
    print(distrb)

    # plt.plot(power_law_dd, label=f'{pgf_formalism.z1_of(power_law_dd)}')
    # plt.plot(power_law_dd10, label=f'{pgf_formalism.z1_of(power_law_dd10)}')
    # plt.plot(power_law_dd100, label=f'{pgf_formalism.z1_of(power_law_dd100)}')
    # plt.plot(poisson, label=f'{pgf_formalism.z1_of(poisson)}')
    # plt.semilogy()
    # plt.xlim([0, 20])
    # plt.ylim([.000001, 1])
    # plt.legend(loc='upper right')
    # plt.show()
    power_law_dd = degree_distributions.power_law_degree_distrb(400)

    # ## Comparative Figures for the talk, with universal intervention

    # plotting_util.phaseSpace_from_data('../data/talk_results/powerlawrollout4565.txt', 5, '$(s,m)$ phase space')
    # plotting_util.phaseSpace_from_data('../data/talk_results/powerlawuniversal45_intv.txt', 5, '$(s,m)$ phase space')

    # plotting_util.plots_for_nerccs_talk([3, 6, 10, 15], 400,  '../data/talk_results/pl_T8_no_intervene.txt', '../data/talk_results/powerlawrollout456{0}.txt',
    #                                                 '../data/talk_results/universal_pl_T8_gen4_T4_intervene.txt', '../data/talk_results/powerlawuniversal4_T4{0}_intv.txt',
    #                                                 same_plot=True, normalize_axis_x=False, plot_distribution_inset=True, inset_label=f'$T_g \\rightarrow T_4=0.4$')
    #
    # ## Universal Intervention
    # # plotting_util.plots_for_nerccs_talk([3, 6, 10, 15], 400,  '../data/talk_results/pl_T8_no_intervene.txt', '../data/talk_results/powerlawrollout456{0}.txt',
    # #                                                 '../data/talk_results/universal_pl_T8_gen4_T4_intervene.txt', '../data/talk_results/powerlawuniversal4_T4{0}_intv.txt',
    # #                                                 same_plot=True, normalize_axis_x=False, inset_label=f'$T_g \\rightarrow T_4=0.4$')
    # # plt.show()
    #
    # ## Random Vaccination
    # plotting_util.plots_for_nerccs_talk([3, 5, 6, 15], 400,  '../data/talk_results/pl_T8_no_intervene.txt', '../data/talk_results/powerlawrollout456{0}.txt',
    #                                                 '../data/talk_results/random_pl_T8_gen4_50pct_intervene.txt', '../data/talk_results/powerlawrandom4{0}_intv.txt',
    #                                                 same_plot=True, normalize_axis_x=False, inset_label=f'$T_g \\rightarrow T_4=0$ for \n 50\% Vaccinated')
    # plt.show()
    #
    # ## Random Rollout vaccination
    # plotting_util.plots_for_nerccs_talk([3, 10, 15, 20], 400,  '../data/talk_results/pl_T8_no_intervene.txt', '../data/talk_results/powerlawrollout456{0}.txt',
    #                                                 '../data/talk_results/rollout_pl_T8_gen456_50pct_intervene.txt', '../data/talk_results/powerlawrollout456{0}_intv.txt',
    #                                                 same_plot=True, normalize_axis_x=False, inset_label=f'$T_g \\rightarrow T_4=0$, for 10\% Vaccinated\n $\\rightarrow T_5=0$, for 30\% Vaccinated\n$\\rightarrow T_6=0$, for 50\% Vaccinated')
    # plt.show()

    ## Theoretical Data for the talk
    # pgf_formalism.compute_phase_space(21, 400, power_law_dd, 0.8, True, [2, 3, 4, 5, 6, 10, 12, 15, 20], '../data/talk_results/powerlawuniversal4_T4{0}',
    #                          intervention_gen=4, intervention_trans=0.4, do_non_interv=False, do_interv=True, intervention_type="universal_intervention")
    # plotting_util.outbreak_size_curves([2, 4, 6, 10, 15], 200, '../data/talk_results/powerlawrollout456{0}.txt', '../data/talk_results/powerlawuniversal4_T4{0}_intv.txt', same_plot=True)

    ### Universal intervention:
    # ensemble.run_ensemble_intervention_effects(power_law_dd, '../data/talk_results/pl_T8', 75000, 10000,
    #                                            init_T=0.8, gen_intervene=4, T_intervene=0.4,
    #                                         intervention_type="none", run_regular=True)
    # plt.figure(figsize=(12,4))
    # #todo axis labels
    # plotting_util.graph_infection_size_distribution_by_gen([2, 4, 6, 10, 15], 200, '../data/talk_results/pl_T8_no_intervene.txt', None
    #                                                   ,'../data/talk_results/universal_pl_T8_gen4_T4_intervene.txt')
    # plt.tight_layout()
    # plt.show()


    # This works and is awesome, so we can just specify which time gens we wanna show:
    # make a plot that has both time and generation side by side, (also save the results) to show the difference in probability distribution
    # it also kind of shows whether generation is a good proxy for time anyways
    # t_buckets, t_distribution = ensemble.ensemble_time_distributions(power_law_dd, 20000, 10000, initial_T=0.6, gamma=0.001)
    #TODO run this tonight, 50k
    # file_root = 'plaw_T8_10k_gamma1_g_over_bg'
    # file_root = 'poiss_T5_5h_q_2_gamma1_g_over_b'
    # # file_root = 'poiss_T5_1k_q_2_gamma1_g_over_b'
    # # file_root = 'poiss_T5_1k_q_2_gamma1_g_over_b'
    # file_root = 'chain_T9_1k_q_1_gamma1_g_over_b'
    # file_root = 'poiss_T8_10k_q_1_gamma1_g_over_b'
    # # ensemble.run_ensemble_with_time_distribution(poisson, f'../data/testing/{file_root}', 30, 10000, 0.8, recover_rate=.001)
    # # todo change to 10k suffix when plotting for wednesday
    # plt.figure('Generations')
    # plt.title('Generational time distribution of cumulative infections')
    # gen_emergence = np.loadtxt(f'../data/testing/{file_root}_gen_emergence_times.txt', delimiter=',')
    # # plotting_util.graph_infection_size_distribution_by_gen([2, 3, 4, 5, 6], 200, f'../data/testing/{file_root}_generational.txt', gen_emergence_times=gen_emergence)
    # # vlines doesn't really make sense here since the x axis is node quantity, but, can use as labels
    # # plt.vlines(gen_emergence, ymin=0, ymax=1, colors='orange', alpha=0.5, linestyles=':')
    # plt.tight_layout()
    # plt.savefig(f'../data/{file_root}_gens.png')
    # plt.figure('Time')
    # t_buckets = np.loadtxt(f'../data/testing/{file_root}_time_distribution_values.txt', delimiter=',')
    # t_distribution = np.loadtxt(f'../data/testing/{file_root}_time_distribution.txt', delimiter=',')
    # colors = ['red', 'orange', 'green', 'blue', 'goldenrod', 'purple', 'pink', 'teal', 'black', 'gold', 'chocolate',
    #           'dodgerblue', 'darkslategray', 'mediumorchid', 'magenta', 'tomato', 'midnightblue',
    #           'cadetblue', 'crimson']
    # for x in [2, 3, 4, 5, 6]:
    #     # beta = .004 #(with T=0.8 and gamma=.001)
    #     beta = .001 #(with T=0.5 and gamma=.001)
    #     # beta = .008 #(with T=.8 and gamma=.002)
    #     plt.plot(t_distribution[x][:2000], label=f't={int(t_buckets[x])}, Expected gen (txBeta) = {int(x)}', color=colors[0], alpha=0.7)
    #     colors.remove(colors[0])
    # plt.legend(loc='upper right')
    # plt.semilogy()
    # plt.xlabel('Cumulative number nodes infected')
    # plt.ylabel('Probability')
    # plt.title('Clock time distribution of cumulative infections')
    # plt.tight_layout()
    # plt.savefig(f'../data/{file_root}_time.png')
    # plt.figure("Generational emergence vs expectation")
    # plt.title('Typical clock time of generation emergence vs expected clock time')
    # max_len = min(len(t_buckets), len(gen_emergence))
    # # max_len = 60
    # max_len = 21
    # plt.scatter(t_buckets[:max_len], gen_emergence[:max_len])
    # plt.xlabel('Expected time of emergence')
    # plt.ylabel('Actual time of emergence')
    # lin_reg_result = stats.linregress(t_buckets[:max_len], gen_emergence[:max_len])
    # x_vals = np.arange(0, 15000)
    # x_vals = np.arange(0, 2000)
    # y_fit = lin_reg_result.intercept + lin_reg_result.slope * x_vals
    # plt.plot(x_vals, y_fit, color='red', label=f'Linear fit of $y={np.round(lin_reg_result.slope, 2)}$')
    # plt.legend(loc='lower right')
    # plt.tight_layout()
    # plt.savefig(f'../data/{file_root}_time_fit.png')
    # plt.show()

    # expected_time = np.zeros(100)
    # for g in range(1, 100):
    #     expected_time[g] = np.sum(expected_time[:g]) + 1/((1.65**g)*beta)
    # plt.scatter(expected_time, gen_emergence)
    # plt.show()
    #
    # correlation = plotting_util.gen_vs_time_correlation(np.loadtxt(f'../data/testing/{file_root}_generational.txt', delimiter=','), t_distribution)
    # plt.scatter(np.arange(len(t_distribution))[2:], correlation[2:])
    # # plt.show()
    #
    # expected_infct = np.loadtxt(f'./../data/expected_cum_0.5_m.txt')
    # expected_time = np.zeros(50)
    # for g in range(len(expected_infct)):
    #     expected_time[g] = g/(beta*expected_infct[g])
    # plt.scatter(expected_time, gen_emergence[:50])
    # plt.show()
    # pgf_formalism.expected_num_infected(power_law_dd, .5)
    # plotting_util.outbreak_size_curves([2, 6, 10, 20], 200, '../../sample/phaseSpaceT6toT4Gen{0}.txt', '../../sample/phaseSpaceT6toT4Gen{0}_intv.txt', same_plot=True)




    # pgf_formalism.compute_phase_space(16, 400, power_law_dd, 0.6, True, [2, 3, 4, 5, 6, 10, 12, 15, 20], '../data/powerlawrollout345{0}',
    #                          4, 0.0, 0.06, rollout_dict={3:0.05, 4:0.10, 5:0.10}, do_non_interv=True, do_interv=True, intervention_type="random_rollout")
    # plotting_util.outbreak_size_curves([4, 6, 10, 15], 400, '../data/randompoisson{0}.txt', '../data/randompoisson{0}_intv.txt', same_plot=True)

    # for g in [3, 4, 5, 6]:
        # plotting_util.phaseSpace_from_data('../data/randomrollout{0}.txt'.format(g), g, 'Power law DD phase space with T=0.6')
        # plotting_util.phaseSpace_from_data('../data/randomrollout{0}_intv.txt'.format(g), g, f'Random Rollout phase space get {g}')

    # Run an ensemble simulation

    # ensemble.run_ensemble_intervention_effects(poisson, '../data/testing/poiss_for_paper', 100, 10000,
    #                                            init_T=0.8, run_regular=True)
    # plotting_util.graph_infection_size_distribution_by_gen([2, 4, 6, 10], 110, '../data/testing/poiss_for_paper_no_intervene.txt'
    #                                                   )
    # # plt.show()
    # plotting_util.plot_sims_vs_analytical_multigens([3, 6, 10], 200,  '../data/testing/powerlawrollout345_no_intervene.txt', '../data/powerlawrollout345{0}.txt',
    #                                                 '../data/testing/powerlawrollout345_intervene.txt', '../data/powerlawrollout345{0}_intv.txt',
    #                                                 same_plot=True, normalize_axis_x=False)
    # # plt.show()


    print('Done')

    # Sample usage of plotting tools
    #
    # plotting_util.plot_sims_vs_analytical_multigens([2, 3, 4], 100, '../../testing/random_vac_40perc_to0gen4_size_distrb_per_gen_no_interv.txt',
    #                                         '../../testing/randomvacT6Gen{0}.txt',
    #                                            '../../testing/random_vac_40perc_to0gen4_size_distrb_per_gen_with_interv.txt',
    #                                            '../../testing/randomvacT6Gen{0}_intv.txt',
    #                                                 same_plot=True, normalize_axis_x=False, plot_distribution_inset=False, grayscale=False)
    # for g in [4, 6]:
    #     plotting_util.phaseSpace_from_data('../../sample/phaseSpaceT6toT4Gen{0}.txt'.format(g), g, 'Power law DD phase space with T=0.6')
    #     plotting_util.phaseSpace_from_data('../../sample/phaseSpaceT6toT4Gen{0}_intv.txt'.format(g), g, 'Power law DD phase space with T=0.6 to T=0.4 at gen 4')
    #
    #
    #
    # plotting_util.plot_sims_vs_analytical_multigens([2, 3, 4], 100, '../../testing/random_vac_40perc_to0gen4_size_distrb_per_gen_no_interv.txt',
    #                                         '../../testing/randomvacT6Gen{0}.txt',
    #                                            '../../testing/random_vac_40perc_to0gen4_size_distrb_per_gen_with_interv.txt',
    #                                            '../../testing/randomvacT6Gen{0}_intv.txt',
    #                                                 same_plot=True, normalize_axis_x=False, plot_distribution_inset=False, grayscale=False)
    #
    # power_law_dd = degree_distributions.power_law_degree_distrb(400)
    # pgf_formalism.phaseSpace(25, 400, power_law_dd, 0.6, True, [2, 3, 4, 5, 6, 10, 15, 20], '../../testing/randomvacT6Gen{0}',
    #                          4, 0.0, 0.4, do_non_interv=False, do_interv=True)

    # plotting_util.graph_infection_size_distribution_by_gen([2, 4, 6, 10], 110, '../../testing/random_vac_40perc_to0gen4_size_distrb_per_gen_no_interv.txt'
    #                                                        ,'../../testing/random_vac_40perc_to0gen4_size_distrb_per_gen_with_interv.txt')
    # plt.show()




    # plotting_util.graph_infection_size_distribution_by_gen([2, 4, 6, 10], 110, '../../testing/random_vac_40perc_gen4_size_distrb_per_gen_no_interv.txt')
    #                                                   #,'../../testing/', 'random_vac_40perc_gen4_size_distrb_per_gen_with_interv.txt')
    #
    #
    # plotting_util.graph_infection_size_distribution_by_gen([2, 4, 6, 10], 110, '../../testing/4fig_power_law_06_to04_gen4_size_distrb_per_gen_no_interv.txt')
    #                                                   #,'../../testing/', '4fig_power_law_06_to04_gen4_size_distrb_per_gen_with_interv.txt')
    # # plt.show()



    # Some more sample usage of plotting functions from data:

    # plotting_util.plot_sims_vs_analytical_multigens([4, 6, 15], 230, '../../testing/4fig_power_law_06_to04_gen4_size_distrb_per_gen_no_interv.txt',
    #                                         '../../sample/phaseSpaceT6toT4Gen{0}.txt',
    #                                            '../../testing/4fig_power_law_06_to04_gen4_size_distrb_per_gen_with_interv.txt',
    #                                            '../../sample/phaseSpaceT6toT4Gen{0}_intv.txt',
    #                                                 same_plot=True, normalize_axis_x=True, plot_distribution_inset=False, grayscale=False)
    #
    #
    # plotting_util.plot_sims_vs_analytical_multigens([2, 4, 10], 230, '../../sample/4fig_power_law_06_to04_gen4_size_distrb_per_gen_no_interv.txt',
    #                                         '../../sample/phaseSpaceT6toT4Gen{0}.txt',
    #                                            '../../sample/4fig_power_law_06_to04_gen4_size_distrb_per_gen_with_interv.txt',
    #                                            '../../sample/phaseSpaceT6toT4Gen{0}_intv.txt',
    #                                                 same_plot=True, normalize_axis_x=True, plot_distribution_inset=True, grayscale=True)

    # # plotting_util.distribution_heatmap(100, 400, pgf_formalism.power_law_degree_distrb(400), 0.5)
    # plotting_util.outbreak_size_curves([2, 6, 10, 20], 200, '../../sample/phaseSpaceT6toT4Gen{0}.txt', '../../sample/phaseSpaceT6toT4Gen{0}_intv.txt', same_plot=True)
    # plotting_util.outbreak_size_curves([2, 6, 10, 20], 200, '../../sample/phaseSpaceT6toT4Gen{0}.txt', None, same_plot=True)
    # plotting_util.outbreak_size_curves([2, 6, 10, 20], 200, '../../sample/phaseSpaceT6toT4Gen{0}.txt', '../../sample/phaseSpaceT6toT4Gen{0}_intv.txt', same_plot=False)
    # plotting_util.outbreak_size_curves([2, 6, 10, 20], 200, '../../sample/phaseSpaceT6toT4Gen{0}.txt', None, same_plot=False)