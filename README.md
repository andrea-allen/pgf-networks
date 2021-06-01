# pgf-networks

This project was created at the University of Vermont, Complex Systems and Data Science Center in the Laboratory for Structure and Dynamics to accompany a working paper.
Contact the authors for specific questions.


# Files included in this repo

## src directory:
### authors' main workspaces:
` main.py`, `main_NJR.py`, `main_AJA.py` are authors' workspaces for running available code. May be used as a reference
but not a consistent source of code or artifacts. Content may change regularly at author's discretion.

### theoretical computation
`degree_distributions.py` - Helper functions for generating certain degree distributions.

`gen_extinct_prob.py` - Code for theoretical extinction probability of an epidemic.

`pgf_formalism.py` - Code for theoretical probability generating functions, and related results.

`plotting_util.py` - Helper functions for plotting results, for comparing empirical results against analytical computations

`figs_for_paper.py` - Specific code for figures for paper draft

`sample_use.py` - Various samples of code snippets for reference

## analysis directory:

`ensemble.py` - Code for running simulations to validate theoretical results. All simulations are run via an original
SIR simulation package `Epintervene` which can be found on GitHub as well. Code in this file `ensemble.py` is code to
aggregate results and run specific validations using the package. For package details, see https://github.com/andrea-allen/epintervene


`covid.py` - Comparing theoretical results to covid data from John's Hopkins.

## test directory:
Directory for unit tests. Currently the unit tests are obsolete/deprecated. Do not use.

## archive directory:
Archive of previous files and deprecated code. Do not use. May be deleted at author's discretion.



The algorithm for the event driven simulations is the Gillespie algorithm.
Gillespie:
- Draw next timestep tau
- Draw the process/event class with rate proportional to the sum total rates of all event classes
    probability of next type of event is proportional to its total rate
- Draw the actual event proportional to the max rate in its rate vector:
using the rejection method:
select an event with uniform chance
accept with probability v_p/v_max, v_max = largest rate in the vector
repeat until an event is accepted


## Theoretical details for pgf_formalism.py

This file contains the analytical tools to perform computations following our mathematical formalism of generating functions.

The collection of methods `pdf_of(), z1_of(), z2_of(), g1_of()` provide the basic computation
for building generating functions and associated properties, given a degree list or degree distribution.

The rest of the file contains algorithms and example code to mainly produce two types of figures:
Phase space of the probability distribution of `s,m` total infected (see our paper for details) and 
the marginal distributions for `s`, the total number infected at each generating `g`.
### Phase Space Computation and Visualization
The method `phaseSpace(num_Gens, num_nodes)` is the logical entrypoint. Users may configure the code to first
produce their desired starting degree distribution. Then this method calls `Psi()` which generates
a 3-dimensional matrix, containing a single 2-d matrix for each generation `g`, containing the probabilities
for each state `(s,m)` as indexed by the rows and columns of the matrix. The method `plot_psi()` will visualize one particular
layer. 

The method `Psi(...)` produces a 3-d set of 2-d matrices, one per generation `g`, with the `(s,m)` phase space distribution.
The methods here follow the recursive algorithm described in the associated paper.
1. We initialize the phase space matrix for generation 0.
1. g_1, g_0 are generated using the given degree distribution, but include
the transmissibility T to create two new generating functions
1. The matrices M_0, M_1 are computed using `constructMatrix(g0, g1)` which generates a matrix in the indicies
`m'` by `m`, where each `m'`th row is the probability distribution for the
convolution of `m'` copies of `g_1` (or `g_0`, appropriately) and the entry in the
`m`th column is the `m`th coefficient of the resulting generating function.
1. Each `s,m`th entry in the current generation's matrix is then computed by
calling `computeLittlePsi` using the recursive algorithm, which requires the appropriate
`M` matrix, the generating functions, the previous generation's `Psi` matrix layer.

### Marginal Outbreak Size Distribution Visualization
The method `outbreak_size_curves(...)` provides the code to compute the marginal probability distribution
for `s` for each generation `g`, which is the probability that there will be `s` infected by generation `g`. 
Calls the `Psi` method, then compresses the rows (which represent `m`) into one and normalizes
the distribution, resulting in the marginal. The bulk of the method is customizable visualiation.

### gen_extinct_prob.py
Code for computing the probability of extinction given an epidemic arrives at some state in generation
`g`, cumulative cases `s`, and current cases `m`. 





