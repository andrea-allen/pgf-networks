from matplotlib import rc
from src import pgf_formalism
from src import SIR_sims
from src import analysis
from src import degree_distributions

if __name__ == '__main__':
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    rc('text', usetex=True)
    print('pgfs yay!')

    # Sample usage:
    # degree_distributions.power_law_degree_distrb()
    power_law_dd = degree_distributions.power_law_degree_distrb(400)
    # Run two sets of ensembles: one base level with no intervention, one with intervention introduced by specified params
    SIR_sims.simulate_intervention_effects(power_law_dd, '../../sample/power_law_06_to04_gen4', 20000, 10000,
                                                               0.6, 4, 0.4, .001)
    pgf_formalism.phaseSpace(25, 400, power_law_dd, 0.8, True, [2, 4, 6, 10, 15, 20], '../../sample/phaseSpaceGen{0}')
    for g in [2, 4, 6, 10, 15, 20]:
        analysis.phaseSpace_from_data('../../sample/phaseSpaceGen{0}.txt'.format(g), g, 'Power law DD phase space with T=0.8')

    # analysis.distribution_heatmap(100, 400, pgf_formalism.power_law_degree_distrb(400), 0.5)
    analysis.outbreak_size_curves([2, 6, 11, 18], 200, '../../pgf-nets-data/allPsiT8_{0}.txt', '../../pgf-nets-data/allPsiT8_{0}_int.txt', same_plot=True)
    analysis.outbreak_size_curves([2, 6, 11, 18], 200, '../../pgf-nets-data/allPsiT8_{0}.txt', None, same_plot=True)
    analysis.outbreak_size_curves([2, 6, 11, 18], 200, '../../pgf-nets-data/allPsiT8_{0}.txt', '../../pgf-nets-data/allPsiT8_{0}_int.txt', same_plot=False)
    analysis.outbreak_size_curves([2, 6, 11, 18], 200, '../../pgf-nets-data/allPsiT8_{0}.txt', None, same_plot=False)
    analysis.plot_sims_vs_analytical_multigens([2, 6, 11, 18], 200, '../../pgf-nets-data/power_law_08_to_06_gen3_size_distrb_per_gen_no_interv.txt',
                                            '../../pgf-nets-data/allPsiT8_{0}.txt',
                                               '../../pgf-nets-data/power_law_08_to_06_gen3_size_distrb_per_gen_with_interv.txt',
                                               '../../pgf-nets-data/allPsiT8_{0}_int.txt',
                                               same_plot=True)

    analysis.plot_sims_vs_analytical_multigens([2, 6, 11, 18], 200, '../../pgf-nets-data/power_law_08_to_06_gen3_size_distrb_per_gen_no_interv.txt',
                                            '../../pgf-nets-data/allPsiT8_{0}.txt',
                                               None,
                                               None,
                                               same_plot=True)

    analysis.plot_sims_vs_analytical_multigens([2, 6, 11, 18], 200,
                                               '../../pgf-nets-data/power_law_08_to_06_gen3_size_distrb_per_gen_no_interv.txt',
                                               '../../pgf-nets-data/allPsiT8_{0}.txt',
                                               None,
                                               None,
                                               same_plot=False)

    analysis.graph_infection_size_distribution_by_gen([2, 3, 4, 6, 15], 250, '../../pgf-nets-data/',
                                                          'power_law_08_to_06_gen3_size_distrb_per_gen_no_interv.txt',
                                                      '../../pgf-nets-data/',
                                                          'power_law_08_to_06_gen3_size_distrb_per_gen_with_interv.txt'
                                                          )

    analysis.graph_infection_size_distribution_by_gen([2, 3, 4, 6, 15], 250, '../../pgf-nets-data/',
                                                          'power_law_08_to04_gen3_fast_size_distrb_per_gen_no_interv.txt',
                                                      '../../pgf-nets-data/',
                                                          'power_law_08_to04_gen3_fast_size_distrb_per_gen_with_interv.txt'
                                                          )
    analysis.phaseSpace_from_data('../../pgf-nets-data/allPsiT8_2.txt', 2, 'Gen 2')
    analysis.phaseSpace_from_data('../../pgf-nets-data/allPsiT8_2_int.txt', 2, 'Gen 2')
    analysis.phaseSpace_from_data('../../pgf-nets-data/allPsiT8_6.txt', 6, 'Gen 6')
    analysis.phaseSpace_from_data('../../pgf-nets-data/allPsiT8_6_int.txt', 6, 'Gen 6')
    analysis.phaseSpace_from_data('../../pgf-nets-data/allPsiT8_18.txt', 18, 'Gen 18')
    analysis.phaseSpace_from_data('../../pgf-nets-data/allPsiT8_18_int.txt', 18, 'Gen 18')


    # TODO add local tests for pgf_formalism to compute correct values, and ensure the saving and reading locally works

