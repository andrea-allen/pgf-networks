from matplotlib import rc
import SIR_sims
import pgf_formalism

if __name__ == '__main__':
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    rc('text', usetex=True)
    print('pgfs yay!')
    # probMatrix = pgf_formalism.phaseSpace(20, 400)
    # pgf_formalism.outbreak_size_curves(20, 400)
    pgf_formalism.plot_sims_vs_analytical_outbreak_sizes()
    SIR_sims.run()
