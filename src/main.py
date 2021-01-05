from matplotlib import rc
from src import pgf_formalism
from src import sandbox

if __name__ == '__main__':
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    rc('text', usetex=True)
    print('pgfs yay!')

    sandbox.run_sandbox()
    # pgf_formalism.newFigure()
    # probMatrix = pgf_formalism.phaseSpace(20, 400)
    # pgf_formalism.outbreak_size_curves(20, 400)
    # pgf_formalism.plot_sims_vs_analytical_outbreak_sizes()
    # SIR_sims.run()

    # pgf_formalism.phaseSpace_from_data('../pgf-nets-data/allPsiT8_2.txt', 2, 'Power Law Degree Distribution')
    # pgf_formalism.phaseSpace_from_data('../pgf-nets-data/allPsiT8_6.txt', 6, 'Power Law Degree Distribution')
    # pgf_formalism.phaseSpace_from_data('../pgf-nets-data/allPsiT8_11.txt', 11, 'Power Law Degree Distribution')
    # pgf_formalism.phaseSpace_from_data('../pgf-nets-data/allPsiT8_18.txt', 18, 'Power Law Degree Distribution')
    #
    # pgf_formalism.phaseSpace_from_data('../pgf-nets-data/allPsiT8_2_int.txt', 2, 'Power Law Degree Distribution')
    # pgf_formalism.phaseSpace_from_data('../pgf-nets-data/allPsiT8_6_int.txt', 6, 'Power Law Degree Distribution')
    # pgf_formalism.phaseSpace_from_data('../pgf-nets-data/allPsiT8_11_int.txt', 11, 'Power Law Degree Distribution')
    # pgf_formalism.phaseSpace_from_data('../pgf-nets-data/allPsiT8_18_int.txt', 18, 'Power Law Degree Distribution')

