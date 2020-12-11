# pgf-networks

# Files included in this repo

## event_driven.py
This file contains the necessary classes and methods to run a single, event-driven,
SIR simulation, with customizable parameters. 

### Entry point:
We recommend understanding `event_driven.py` functionality by observing the class
`Simulation`, and its arguments in its constructor. Once a Simulation object is created,
one need only call `Simulation.run_sim()` with its respective arguments to run one realization
of the simulation. Results are then stored in the object itself and can be accessed by calling
`Simulation.generate_matrix_gen(num_gens)` to obtain a 2-row list of `s` and `m` per generation,
the number `s` infected overall during that generation and the number `m` infected only during the most
recent generation. 

Suggested use:
TODO rename matrix method
```
G = #must be a networkX graph
sim = event_driven.Simulation(1000000, G, Lambda, Gamma, pos)
    sim.run_sim(intervention_gen, beta_interv)
    sm_matrix = sim.generate_matrix_gen(20)
```

There are many more intricacies in the `event_driven.py` file, which will be documented 
at a future date. Most of them can be inferred from their usage in the code.

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

## main.py
Driver for calling methods from the repository to run simulations, read back
data, produce figures, etc. 



Given: Network adjacency matrix
Outcome: Number infected over time
Book-keeping: Generations of infection, based on nodes tagged

data structure:
gens_infected: [1, 3, 7, 9, ...]

V_IS = edge list of infected node id's to susceptible node ID's
V_IS_p = positions of the susceptible vertices
L_p = rates of each event happening (lambda_ij)
total infection, compute after each event
L = sum_ij^N lambda_ij delta(signma_i, 1)delta(sigma_j, 0) = sum_p=1^N_IS L_p
V_I = id's of infected nodes, V_I_p position of infected nodes
M_p = recovery rates
delta_i=0 S, =1 I, =2 R
When an event happens and a node changes class, move the pth node out, move the last node to position p and shorten the list
V_R = recovered nodes
A = waning immunity rates

Total rate of recovery:
M = sum_i=1^N mu_i delta(sigma_i, 1) = sum_p=1^N_I M_p
Total rate of waning immunity:
A = sum_i=1^N alpha_i delta(sigma_i, 2) = sum_p=1^N_R A_p
At each event, update M or A by adding or subtracting from M_p or A_p when new nodes are added or removed to V_I, V_R


Gillespie:
- Draw next timestep tau
- Draw the process/event class with R = sum total rates of all event classes
    probability of next type of event is proportional to its total rate
- Draw the actual event proportional to the max rate in its rate vector:
using the rejection method:
select an event with uniform chance
accept with probability v_p/v_max, v_max = largest rate in the vector
repeat until an event is accepted
for simplicity, let the max's be global over the whole network

Our twist:
include a vector of V_I_g which includes its generation
when an event happens from V_IS, add S to V_I, and add V_I_g +1 for that new node based on I
would be good if the nodes were object oriented
