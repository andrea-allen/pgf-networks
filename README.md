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
