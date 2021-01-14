from matplotlib import rc
from src import pgf_formalism
from src import SIR_sims
from src import analysis

if __name__ == '__main__':
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    rc('text', usetex=True)
    print('pgfs yay!')

    # Sample usage:
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


    # pgf_formalism.newFigure()
    # probMatrix = pgf_formalism.phaseSpace(20, 400)
    # SIR_sims.run()

