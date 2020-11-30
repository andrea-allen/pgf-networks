from matplotlib import rc
import SIR_sims
import pgf_formalism

if __name__ == '__main__':
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    rc('text', usetex=True)
    print('pgfs yay!')
    probMatrix = pgf_formalism.phaseSpace(20)
    SIR_sims.run()
