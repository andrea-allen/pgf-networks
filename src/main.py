import time

import matplotlib.pyplot as plt
from matplotlib import rc

from analysis import ensemble
from src import degree_distributions
from src import pgf_formalism
from src import plotting_util

if __name__ == '__main__':
    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)
    print('pgfs yay!')

    # Run an ensemble simulation
    power_law_dd = degree_distributions.power_law_degree_distrb(400)

    ensemble.simulate_intervention_effects(power_law_dd, '../data/testing/multi_intervention_tiny', 20000, 2000,
                                                               0.6, 4, 0.0, .001)
    plotting_util.graph_infection_size_distribution_by_gen([2, 4, 6, 10], 110, '../data/testing/multi_intervention_tiny_no_intervene.txt'
                                                      ,'../data/testing/multi_intervention_tiny_intervene.txt')
    plt.show()


    print('Done')

    # Sample usage of plotting tools

    # plotting_util.plot_sims_vs_analytical_multigens([2, 3, 4], 100, '../../testing/random_vac_40perc_to0gen4_size_distrb_per_gen_no_interv.txt',
    #                                         '../../testing/randomvacT6Gen{0}.txt',
    #                                            '../../testing/random_vac_40perc_to0gen4_size_distrb_per_gen_with_interv.txt',
    #                                            '../../testing/randomvacT6Gen{0}_intv.txt',
    #                                                 same_plot=True, normalize_axis_x=False, plot_distribution_inset=False, grayscale=False)
    # for g in [4, 6]:
    #     analysis.phaseSpace_from_data('../../sample/phaseSpaceT6toT4Gen{0}.txt'.format(g), g, 'Power law DD phase space with T=0.6')
    #     analysis.phaseSpace_from_data('../../sample/phaseSpaceT6toT4Gen{0}_intv.txt'.format(g), g, 'Power law DD phase space with T=0.6 to T=0.4 at gen 4')



    plotting_util.plot_sims_vs_analytical_multigens([2, 3, 4], 100, '../../testing/random_vac_40perc_to0gen4_size_distrb_per_gen_no_interv.txt',
                                            '../../testing/randomvacT6Gen{0}.txt',
                                               '../../testing/random_vac_40perc_to0gen4_size_distrb_per_gen_with_interv.txt',
                                               '../../testing/randomvacT6Gen{0}_intv.txt',
                                                    same_plot=True, normalize_axis_x=False, plot_distribution_inset=False, grayscale=False)

    power_law_dd = degree_distributions.power_law_degree_distrb(400)
    # pgf_formalism.phaseSpace(25, 400, power_law_dd, 0.6, True, [2, 3, 4, 5, 6, 10, 15, 20], '../../testing/randomvacT6Gen{0}',
    #                          4, 0.0, 0.4, do_non_interv=False, do_interv=True)
    mytime = time.time()
    # SIR_sims.simulate_intervention_effects(power_law_dd, '../../testing/random_vac_40perc_to0gen4', 20000, 10000,
    #                                                            0.6, 4, 0.0, .001)
    print('Total, ', time.time()-mytime)
    plotting_util.graph_infection_size_distribution_by_gen([2, 4, 6, 10], 110, '../../testing/random_vac_40perc_to0gen4_size_distrb_per_gen_no_interv.txt'
                                                           ,'../../testing/random_vac_40perc_to0gen4_size_distrb_per_gen_with_interv.txt')
    plt.show()




    plotting_util.graph_infection_size_distribution_by_gen([2, 4, 6, 10], 110, '../../testing/random_vac_40perc_gen4_size_distrb_per_gen_no_interv.txt')
                                                      #,'../../testing/', 'random_vac_40perc_gen4_size_distrb_per_gen_with_interv.txt')


    plotting_util.graph_infection_size_distribution_by_gen([2, 4, 6, 10], 110, '../../testing/4fig_power_law_06_to04_gen4_size_distrb_per_gen_no_interv.txt')
                                                      #,'../../testing/', '4fig_power_law_06_to04_gen4_size_distrb_per_gen_with_interv.txt')
    plt.show()



    # Some more sample usage of plotting functions from data:

    plotting_util.plot_sims_vs_analytical_multigens([4, 6, 15], 230, '../../testing/4fig_power_law_06_to04_gen4_size_distrb_per_gen_no_interv.txt',
                                            '../../sample/phaseSpaceT6toT4Gen{0}.txt',
                                               '../../testing/4fig_power_law_06_to04_gen4_size_distrb_per_gen_with_interv.txt',
                                               '../../sample/phaseSpaceT6toT4Gen{0}_intv.txt',
                                                    same_plot=True, normalize_axis_x=True, plot_distribution_inset=False, grayscale=False)


    plotting_util.plot_sims_vs_analytical_multigens([2, 4, 10], 230, '../../sample/4fig_power_law_06_to04_gen4_size_distrb_per_gen_no_interv.txt',
                                            '../../sample/phaseSpaceT6toT4Gen{0}.txt',
                                               '../../sample/4fig_power_law_06_to04_gen4_size_distrb_per_gen_with_interv.txt',
                                               '../../sample/phaseSpaceT6toT4Gen{0}_intv.txt',
                                                    same_plot=True, normalize_axis_x=True, plot_distribution_inset=True, grayscale=True)

    # # analysis.distribution_heatmap(100, 400, pgf_formalism.power_law_degree_distrb(400), 0.5)
    # analysis.outbreak_size_curves([2, 6, 10, 20], 200, '../../sample/phaseSpaceT6toT4Gen{0}.txt', '../../sample/phaseSpaceT6toT4Gen{0}_intv.txt', same_plot=True)
    # analysis.outbreak_size_curves([2, 6, 10, 20], 200, '../../sample/phaseSpaceT6toT4Gen{0}.txt', None, same_plot=True)
    # analysis.outbreak_size_curves([2, 6, 10, 20], 200, '../../sample/phaseSpaceT6toT4Gen{0}.txt', '../../sample/phaseSpaceT6toT4Gen{0}_intv.txt', same_plot=False)
    # analysis.outbreak_size_curves([2, 6, 10, 20], 200, '../../sample/phaseSpaceT6toT4Gen{0}.txt', None, same_plot=False)