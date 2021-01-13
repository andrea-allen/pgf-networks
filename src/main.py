from matplotlib import rc
from src import pgf_formalism
from src import sandbox
from src import SIR_sims
from src import analysis

if __name__ == '__main__':
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    rc('text', usetex=True)
    print('pgfs yay!')

    # Sample usage:
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

    # sandbox.run_sandbox()
    # pgf_formalism.newFigure()
    # probMatrix = pgf_formalism.phaseSpace(20, 400)
    # pgf_formalism.outbreak_size_curves(20, 400)
    # pgf_formalism.plot_sims_vs_analytical_outbreak_sizes()
    SIR_sims.run()

    # pgf_formalism.phaseSpace_from_data('../pgf-nets-data/allPsiT8_2.txt', 2, 'Power Law Degree Distribution')
    # pgf_formalism.phaseSpace_from_data('../pgf-nets-data/allPsiT8_6.txt', 6, 'Power Law Degree Distribution')
    # pgf_formalism.phaseSpace_from_data('../pgf-nets-data/allPsiT8_11.txt', 11, 'Power Law Degree Distribution')
    # pgf_formalism.phaseSpace_from_data('../pgf-nets-data/allPsiT8_18.txt', 18, 'Power Law Degree Distribution')
    #
    # pgf_formalism.phaseSpace_from_data('../pgf-nets-data/allPsiT8_2_int.txt', 2, 'Power Law Degree Distribution')
    # pgf_formalism.phaseSpace_from_data('../pgf-nets-data/allPsiT8_6_int.txt', 6, 'Power Law Degree Distribution')
    # pgf_formalism.phaseSpace_from_data('../pgf-nets-data/allPsiT8_11_int.txt', 11, 'Power Law Degree Distribution')
    # pgf_formalism.phaseSpace_from_data('../pgf-nets-data/allPsiT8_18_int.txt', 18, 'Power Law Degree Distribution')

