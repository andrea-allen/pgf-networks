from matplotlib import rc
from src import pgf_formalism
from src import SIR_sims
from src import analysis
from src import degree_distributions
import numpy as np
import matplotlib.pyplot as plt
from src import temporal_wip
import time

if __name__ == '__main__':
    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)
    print('pgfs yay!')

    # # Run two sets of ensembles: one base level with no intervention, one with intervention introduced by specified params
    power_law_dd = degree_distributions.power_law_degree_distrb()
    mytime = time.time()
    SIR_sims.simulate_intervention_effects(power_law_dd, '../../testing/4fig_power_law_06_to04_gen4', 10000, 1000,
                                                               0.6, 4, 0.4, .001)
    print('Total, ', time.time()-mytime)
    analysis.graph_infection_size_distribution_by_gen([2, 6, 15], 100, '../../testing/', '4fig_power_law_06_to04_gen4_size_distrb_per_gen_no_interv.txt',
                                                      '../../testing/', '4fig_power_law_06_to04_gen4_size_distrb_per_gen_with_interv.txt')
    analysis.plot_sims_vs_analytical_multigens([4, 6, 15], 230, '../../testing/4fig_power_law_06_to04_gen4_size_distrb_per_gen_no_interv.txt',
                                            '../../sample/phaseSpaceT6toT4Gen{0}.txt',
                                               '../../testing/4fig_power_law_06_to04_gen4_size_distrb_per_gen_with_interv.txt',
                                               '../../sample/phaseSpaceT6toT4Gen{0}_intv.txt',
                                               same_plot=True, normalize_axis_x=True, plot_distribution_inset=False, grayscale=False)

    # Playing with interventions:
    temporal_wip.simulate_two_layers()

    # Sample usage:
    power_law_dd = degree_distributions.power_law_degree_distrb(400, 2, 1000)
    # print('Mu 1000, Alpha 2 Threshold', degree_distributions.compute_T_threshold_powerlaw(power_law_dd, 1000))
    plt.plot(power_law_dd[:15])
    plt.semilogy()
    # plt.show()
    print('Mu 1000, Alpha 2 mean', degree_distributions.mean_degree(power_law_dd))

    power_law_dd = degree_distributions.power_law_degree_distrb(400, 2, 5)
    print('Mu 5 Alpha 2 Threshold', degree_distributions.compute_T_threshold_powerlaw(power_law_dd, 5))
    plt.plot(power_law_dd[:15])
    plt.semilogy()

    power_law_dd = degree_distributions.power_law_degree_distrb(400, 1.1, 10)
    print('Mu 10 Alpha 1.1 Mean degree', degree_distributions.mean_degree(power_law_dd))
    print('Mu 10 Alpha 1.1 threshold', degree_distributions.compute_T_threshold_powerlaw(power_law_dd, 10))
    plt.plot(power_law_dd[:15])

    plt.show()
    #Note that power law simulations end up still with not many huge outbreaks even with networks 10x the size, is this expected result
    # TODO might need to have predictions go out further? ask team next week
    print(degree_distributions.mean_degree(power_law_dd))
    binomial_dd = degree_distributions.binomial_degree_distb(400)
    print(degree_distributions.mean_degree(binomial_dd))
    print('z1, should be same as mean', pgf_formalism.z1_of(binomial_dd))
    print('Threshold ', degree_distributions.compute_T_threshold_binom(binomial_dd))

    # # Run two sets of ensembles: one base level with no intervention, one with intervention introduced by specified params
    SIR_sims.simulate_intervention_effects(power_law_dd, '../../testing/4fig_power_law_06_to04_gen4', 1000, 10000,
                                                               0.6, 4, 0.4, .001)

    # Theoretical predictions curves match better with the simulation data on the sims with 10k nodes as opposed
    # to 1k nodes. Probably because of the finite size effects. curves look better than they did in our paper, for example
    analysis.plot_sims_vs_analytical_multigens([4, 6, 15], 230, '../../sample/4fig_power_law_06_to04_gen4_size_distrb_per_gen_no_interv.txt',
                                            '../../sample/phaseSpaceT6toT4Gen{0}.txt',
                                               '../../sample/4fig_power_law_06_to04_gen4_size_distrb_per_gen_with_interv.txt',
                                               '../../sample/phaseSpaceT6toT4Gen{0}_intv.txt',
                                               same_plot=True, normalize_axis_x=True, plot_distribution_inset=False, grayscale=False)


    # pgf_formalism.phaseSpace(25, 400, power_law_dd, 0.6, True, [2, 4, 6, 10, 15, 20], '../../sample/phaseSpaceT6Gen{0}', 4, 0.4)
    # pgf_formalism.phaseSpace(25, 400, power_law_dd, 0.6, True, [2, 4, 6, 10, 15, 20], '../../sample/phaseSpaceT6toT4Gen{0}', 4, 0.4)
    # for g in [2, 4, 6, 10, 15, 20]:
    #     analysis.phaseSpace_from_data('../../sample/phaseSpaceT6toT4Gen{0}.txt'.format(g), g, 'Power law DD phase space with T=0.6')
    #     analysis.phaseSpace_from_data('../../sample/phaseSpaceT6toT4Gen{0}_intv.txt'.format(g), g, 'Power law DD phase space with T=0.6 to T=0.4 at gen 4')
    #
    # # analysis.distribution_heatmap(100, 400, pgf_formalism.power_law_degree_distrb(400), 0.5)
    # analysis.outbreak_size_curves([2, 6, 10, 20], 200, '../../sample/phaseSpaceT6toT4Gen{0}.txt', '../../sample/phaseSpaceT6toT4Gen{0}_intv.txt', same_plot=True)
    # analysis.outbreak_size_curves([2, 6, 10, 20], 200, '../../sample/phaseSpaceT6toT4Gen{0}.txt', None, same_plot=True)
    # analysis.outbreak_size_curves([2, 6, 10, 20], 200, '../../sample/phaseSpaceT6toT4Gen{0}.txt', '../../sample/phaseSpaceT6toT4Gen{0}_intv.txt', same_plot=False)
    # analysis.outbreak_size_curves([2, 6, 10, 20], 200, '../../sample/phaseSpaceT6toT4Gen{0}.txt', None, same_plot=False)
    # analysis.plot_sims_vs_analytical_multigens([2, 6, 10, 20], 200, '../../sample/power_law_06_to04_gen4_size_distrb_per_gen_no_interv.txt',
    #                                         '../../sample/phaseSpaceT6toT4Gen{0}.txt',
    #                                            '../../sample/power_law_06_to04_gen4_size_distrb_per_gen_with_interv.txt',
    #                                            '../../sample/phaseSpaceT6toT4Gen{0}_intv.txt',
    #                                            same_plot=True)
    #
    # analysis.plot_sims_vs_analytical_multigens([2, 6, 11, 18], 200, '../../pgf-nets-data/power_law_08_to_06_gen3_size_distrb_per_gen_no_interv.txt',
    #                                         '../../pgf-nets-data/allPsiT8_{0}.txt',
    #                                            None,
    #                                            None,
    #                                            same_plot=True)
    #
    # analysis.plot_sims_vs_analytical_multigens([2, 6, 11, 18], 200,
    #                                            '../../pgf-nets-data/power_law_08_to_06_gen3_size_distrb_per_gen_no_interv.txt',
    #                                            '../../pgf-nets-data/allPsiT8_{0}.txt',
    #                                            None,
    #                                            None,
    #                                            same_plot=False)
    #
    # analysis.graph_infection_size_distribution_by_gen([2, 3, 4, 6, 15], 250, '../../pgf-nets-data/',
    #                                                       'power_law_08_to_06_gen3_size_distrb_per_gen_no_interv.txt',
    #                                                   '../../pgf-nets-data/',
    #                                                       'power_law_08_to_06_gen3_size_distrb_per_gen_with_interv.txt'
    #                                                       )
    #
    # analysis.graph_infection_size_distribution_by_gen([2, 3, 4, 6, 15], 250, '../../pgf-nets-data/',
    #                                                       'power_law_08_to04_gen3_fast_size_distrb_per_gen_no_interv.txt',
    #                                                   '../../pgf-nets-data/',
    #                                                       'power_law_08_to04_gen3_fast_size_distrb_per_gen_with_interv.txt'
    #                                                       )
    # analysis.phaseSpace_from_data('../../pgf-nets-data/allPsiT8_2.txt', 2, 'Gen 2')
    # analysis.phaseSpace_from_data('../../pgf-nets-data/allPsiT8_2_int.txt', 2, 'Gen 2')
    # analysis.phaseSpace_from_data('../../pgf-nets-data/allPsiT8_6.txt', 6, 'Gen 6')
    # analysis.phaseSpace_from_data('../../pgf-nets-data/allPsiT8_6_int.txt', 6, 'Gen 6')
    # analysis.phaseSpace_from_data('../../pgf-nets-data/allPsiT8_18.txt', 18, 'Gen 18')
    # analysis.phaseSpace_from_data('../../pgf-nets-data/allPsiT8_18_int.txt', 18, 'Gen 18')