# pgf-networks
#TODO edit this next
## Readme not currently up to date, check back soon.


This project was created at the University of Vermont, Complex Systems and Data Science Center in the Laboratory for Structure and Dynamics to accompany a working paper titled
"Disease Interventions in Probabilistic Modeling: Expansion and Modeling of
Time Evolution Probability Generating Function Formalism" and contains the original code
used to create our analytical figures, computational results, and event-driven
simulations. 

We recommend exploring the code through trial and error and experimentation, and referencing the usage in the main codebase. Below is
a brief overview of the components of the codebase. Contact the authors for specific questions.


# Files included in this repo

## main.py
Driver for calling methods from the repository to run simulations, read back
data, produce figures, etc. 

## event_driven.py
This file contains the necessary classes and methods to run a single, event-driven,
SIR simulation, with customizable parameters. 

### Entry point:
We recommend understanding `event_driven.py` functionality by observing the class
`Simulation`, and its arguments in its constructor. Once a Simulation object is created,
one need only call `Simulation.run_sim()` with its respective arguments to run one realization
of the simulation. Results are then stored in the object itself and can be accessed by calling
`Simulation.total_infect_over_all_gens(num_gens)` to obtain a 2-row list of `s` and `m` per generation,
the number `s` infected overall during that generation and the number `m` infected only during the most
recent generation. 

Suggested use:
```
G = #must be a networkX graph
sim = event_driven.Simulation(1000000, G, Lambda, Gamma, pos)
    sim.run_sim(intervention_gen, beta_interv)
    results = sim.total_infect_over_all_gens(20)
```

There are many more intricacies in the `event_driven.py` file, which will be documented 
at a future date. Most of them can be inferred from their usage in the code.

The algorithm for the event driven simulations is the Gillespie algorithm.
Gillespie:
- Draw next timestep tau
- Draw the process/event class with R = sum total rates of all event classes
    probability of next type of event is proportional to its total rate
- Draw the actual event proportional to the max rate in its rate vector:
using the rejection method:
select an event with uniform chance
accept with probability v_p/v_max, v_max = largest rate in the vector
repeat until an event is accepted

### Details

####Class `Node`: 
Represents a single node. Contains attributes for the node's label,
state (infectious, susceptible, or recovered), recovery rate, and generation of infection.

####Class `Edge`: 
Contains two nodes, a left and a right, each consisting of a `Node` class object.
Contains methods to infect the right-node and set its generation based on the generation of the 
left-node. 

####Class `Simulation`: 
Object to house and run all relevant information for a single simulation of the event-driven SIR model.

`V_IS`: List of current infected-susceptible edges, changes throughout simulation

`V_I`: List of infected nodes, changes throughout simulation

`run_sim()`: Runs one simulation, utilizes other associated methods within the `Simulation` class

`visualize_network()`: If `visualize` is set to True in the `run_sim()` method, then
at configurable steps of the simulation, a NetworkX graph object will be generated
and the current state of the network will be visualized.

`total_infect_over_all_gens()`: Should be called after the completion of a simulation run,
will return a two-tier vector containing one time series (where time is in terms of generations of infection)
in terms of `m`, the number infected per generation, and `s`, the total number infected by that generation.

`intervene()`: Method that is utilized if the Simulation is initialized with an intervention generation and 
transmission probability change

`draw_tau(), draw_event_class(), draw_specific_event(), determine_draw_tau()`: Helper methods
involved in stochastically yielding the next event time and type.




## SIR_sims.py
This is the driver for running ensembles of simulations, saving and plotting the results,
producing figures, etc. Most of the code in this file is configurable for the user's
desired results, making it difficult to annotate without explicitly pointing out where
certain parameters and filenames are specific to the original uses. 

The method `run()` is just a driver for whatever internal methods the user wishes to call to run and save or plot simulation data.

The method `generate_graph(N, deg_dist)` will return a NetworkX object, G, and 
a drawing position `pos` in order that the same position be saved for drawing
purposes during visualization in the evnt-driven simulations. 
Input the `generate_graph` is `N`, the number of nodes, and the degree 
distribution as a proper probability distribution of any length. Graphs
will then be randomly drawn using the configuration model and the provided
degree distribution. 

The method `simulate(G, pos, Lambda, Gamma, current, intervention_gen, beta_interv)`
is the driving method to return a two-tier vector of results, where there are
`g` generation columns and 2 rows per column, encoding the `m` number of infected
nodes during that generation, and the `s` total number infected at that generation
in the simulation. `Lambda` is an `N x N` matrix encoding the infection probability `beta_i,j` for 
every pair of nodes `i,j`, and can be configured any way the user wants before passing it in.
`Gamma` similarly is a length `N` list of the indiviudal recovery rates for every node.
In the classic SIR models, we keep `Lambda` and `Gamma` values completely
regular for all nodes and node pairs. However, the configurability makes it possible
for the user to classify different nodes and node pairs by different infection
and recovery rates. 

The methods `outbreak_size_distrb_per_gen(...)` and `outbreak_size_distrb_per_gen_with_intervention(...)`
are by and large the same, they run an ensemble of simulations using the
`simulate()` method and the given configurations and then return
a probability distribution, for each generation, of the total outbreak size (`s`)
at that generation. The format of the returned matrix is in `g` rows, each row is a distribution
over `s`, the values of the columns. Arguments required by this method are the
degree distribution, the number of simulations for the ensemble, the number of nodes `N`, the transmission
probability over an edge (`T`) and recovery probability `gamma` from which `beta` will be inferred automatically.

The method `read_back_data()` is another configurable method purely for the user
to organize what data they want to read back after it has been saved for the purposes
of visualization. Some example usage is shown.

The method `simulate_and_compare_rounds_with_with_without_intervention()` similarly is a placeholder
that is configurable, used for our original purposes to run specific simulations.





## pgf_formalism.py

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


### Comparing Empirical Data With Analytical Results
The method `plot_sims_vs_analytical_outbreak_sizes()` is also configurable by the user, the function of this method
is to plot the analytical marginal distribution curves `p_s_g` for particular generations
alongside the empirical results (documented above) from `SIR_sims` that have been previously saved into some data file.

It is worth mentioning that several points in the current version of the code either
read back data for the phase space and marginals, but can be configured to re-compute them from scratch
for a user without saved files. 


