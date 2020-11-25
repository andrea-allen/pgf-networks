# pgf-networks


Andrea's addition to the README demo


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
