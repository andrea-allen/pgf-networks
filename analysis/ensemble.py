import networkx as nx
import numpy as np
from epintervene.simobjects import simulation
from epintervene.simobjects import extended_simulation
from epintervene.simobjects import network
import time


def run_ensemble_intervention_effects(degree_distrb, base_file_name='sim_results',
                                      num_sims=10000, num_nodes=1000, init_T=0.8,
                                      gen_intervene=3, T_intervene=0.4, recover_rate=.001, prop_reduced=0.0,
                                      intervention_gen_list=None, beta_redux_list=None, prop_reduced_list=None,
                                      intervention_type="none",
                                      run_regular=True):
    # Runs both regular and intervention ensembles with the same starting parameters

    # Ensemble run with no intervention:
    if run_regular:
        start_time = time.time()
        size_distrb_per_gen_no_intervention = simulate_ensemble(degree_distrb, num_sims, num_nodes,
                                                                -1, 0.0, init_T, recover_rate, 0.0,
                                                                intervention_gen_list=None, beta_redux_list=None,
                                                                prop_reduced_list=None,
                                                                intervention_type="none"
                                                                )
        print('Total time for ensemble was ', time.time() - start_time, ' with N=', num_nodes, ' nodes and ', num_sims,
              ' number of simulations.')
        np.savetxt(base_file_name + '_no_intervene.txt',
                   size_distrb_per_gen_no_intervention,
                   delimiter=',')

    # Ensemble run with intervention as specified:
    start_time = time.time()
    size_distrb_per_gen_intervention = simulate_ensemble(degree_distrb, num_sims, num_nodes,
                                                         gen_intervene, T_intervene, init_T,
                                                         recover_rate, prop_reduced,
                                                         intervention_gen_list, beta_redux_list,
                                                         prop_reduced_list,
                                                         intervention_type
                                                         )
    print('Total time for ensemble was ', time.time() - start_time, ' with N=', num_nodes, ' nodes and ', num_sims,
          ' number of simulations.')
    np.savetxt(base_file_name + '_intervene.txt', size_distrb_per_gen_intervention,
               delimiter=',')


# Simulate ensemble saves results for generational time series, as opposed to real time results. Those are TBD
def simulate_ensemble(degree_distrb, num_sims=10, N=1000, intervention_gen=-1, intervention_T=0.0, initial_T=0.8,
                      gamma=0.1, prop_reduced=0.0, intervention_gen_list=None, beta_redux_list=None, prop_reduced_list=None,
             intervention_type="none"):
    # Configuring the parameters
    beta_init = -(gamma * initial_T) / (initial_T - 1)
    beta_interv = -(gamma * intervention_T) / (intervention_T - 1)

    # Setting up a results data structure
    outbreak_size_distrb_per_gen_matrix = np.zeros((100, N))

    # Generating the initial network
    G, pos = network.NetworkBuilder.from_degree_distribution(N, degree_distrb)
    A = np.array(nx.adjacency_matrix(G).todense())
    num_nodes_in_net = len(A[0])

    # Generating the infection and recovery rate values
    Beta = np.full((num_nodes_in_net, num_nodes_in_net), beta_init)
    Gamma = np.full(num_nodes_in_net, gamma)

    # Running the ensemble
    for i in range(num_sims):
        # Generate a new network every 500 simulations
        if i % 500 == 0:
            G, pos = network.NetworkBuilder.from_degree_distribution(N, degree_distrb)
            A = np.array(nx.adjacency_matrix(G).todense())
            num_nodes_in_net = len(A[0])
            Beta = np.full((num_nodes_in_net, num_nodes_in_net), beta_init)
            Gamma = np.full(num_nodes_in_net, gamma)
        # Get results of type generational time series
        results = run_single_simulation(A, Beta, Gamma, i, 'generation', intervention_gen, beta_interv, prop_reduced,
                                        intervention_gen_list, beta_redux_list, prop_reduced_list, intervention_type)
        # Recording ensemble results
        for gen in range(len(results)):
            num_total_infctd = int(results[gen])
            try:
                outbreak_size_distrb_per_gen_matrix[gen][num_total_infctd] += 1
            except IndexError:
                print('Index error for g: ', gen, ', gen_s: ', num_total_infctd)
                continue
    # averaging results:
    for gen in range(100):
        gen_time_series = outbreak_size_distrb_per_gen_matrix[gen]
        gen_time_series = gen_time_series / num_sims
        outbreak_size_distrb_per_gen_matrix[gen] = gen_time_series
    return outbreak_size_distrb_per_gen_matrix


def run_single_simulation(A, Beta, Gamma, current, results_type='generation', intervention_gen=-1, beta_interv=-1.0,
                          prop_reduced=0.0, intervention_gen_list=None, beta_redux_list=None, prop_reduced_list=None,
                          intervention_type="none"):
    start_time = time.time()
    # Constructing the simulation of specified Intervention Type, otherwise will run a regular simulation
    if intervention_type == "random-rollout":
        sim = extended_simulation.MultiInterventionSim(A)
        sim.add_infection_event_rates(Beta)
        sim.add_recover_event_rates(Gamma)
        sim.configure_intervention(intervention_gen_list=intervention_gen_list, beta_redux_list=beta_redux_list,
                                   proportion_reduced_list=prop_reduced_list)

    elif intervention_type == "random":
        sim = extended_simulation.RandomInterventionSim(A)
        sim.add_infection_event_rates(Beta)
        sim.add_recover_event_rates(Gamma)
        sim.configure_intervention(intervention_gen=intervention_gen, beta_redux=beta_interv,
                                   proportion_reduced=prop_reduced)
    else:
        sim = simulation.Simulation(A)
        sim.add_infection_event_rates(Beta)
        sim.add_recover_event_rates(Gamma)

    # Run the simulation
    sim.run_sim()

    # Printing progress of the ensemble for the user to keep track
    if current % 50 == 0:
        print('current sim ' + str(current))
        print("--- %s seconds to run simulation---" % (time.time() - start_time))

    # Tabulating results based on user's specified input type
    if str(results_type).lower() == 'generation':
        generational_results = sim.tabulate_generation_results(100)
        return generational_results
    elif str(results_type).lower() == 'time':
        timeseries, timeseries_results_inf, timeseries_results_rec = sim.tabulate_continuous_time(1000)
        return timeseries, timeseries_results_inf, timeseries_results_rec
    elif str(results_type).lower() == 'time_groups':
        timeseries, infection_timeseries_groups = sim.tabulate_continuous_time_with_groups(1000)
        return timeseries, infection_timeseries_groups
    else:
        print('No results type provided, returning basic time series results')
        timeseries, timeseries_results_inf, timeseries_results_rec = sim.tabulate_continuous_time(1000)

    return timeseries, timeseries_results_inf


##### In PROGRESS ######
def ensemble_time_series(network, Beta, Gamma, numsims=100, N=1000):
    A = np.array(nx.adjacency_matrix(network).todense())
    num_nodes_in_net = len(A[0])
    # Beta = np.full((num_nodes_in_net, num_nodes_in_net), beta_init)
    # Gamma = np.full(num_nodes_in_net, gamma)

    total_ts = np.zeros((numsims, 1000))
    total_ts_rec = np.zeros((numsims, 1000))
    ts = np.zeros(1000)

    for i in range(numsims):
        timeseries, timeseries_results_inf, timeseries_results_rec = run_single_simulation(A, Beta, Gamma, i, 'time')
        ts = timeseries
        total_ts[i] = timeseries_results_inf
        total_ts_rec[i] = timeseries_results_rec

    infection_ts = np.mean(total_ts, axis=0)
    recover_ts = np.mean(total_ts_rec, axis=0)
    return ts, infection_ts, recover_ts


def ensemble_time_series_groups(network, Beta, Gamma, numsims=100, N=1000):
    A = np.array(nx.adjacency_matrix(network).todense())
    num_nodes_in_net = len(A[0])
    # Beta = np.full((num_nodes_in_net, num_nodes_in_net), beta_init)
    # Gamma = np.full(num_nodes_in_net, gamma)
    staff_pop = int(len(A) / 3)
    res_pop = int(len(A) / 3)
    county_pop = int(len(A) / 3)

    node_membership_vector = []
    for i in range(staff_pop):
        node_membership_vector.append('staff')
    for j in range(res_pop):
        node_membership_vector.append('residents')
    for k in range(county_pop):
        node_membership_vector.append('county')
    sim = simulation.Simulation(A, membership_groups=['staff', 'residents', 'county'],
                                node_memberships=node_membership_vector)

    total_ts = np.zeros((numsims, 1000))
    total_ts_rec = np.zeros((numsims, 1000))
    ts = np.zeros(1000)

    for i in range(numsims):
        timeseries, timeseries_results_inf, timeseries_results_rec = run_single_simulation(A, Beta, Gamma, i, 'time_groups')
        ts = timeseries  # TODO is timeseries always the same? no, because it splits the max time into bins so maybe still not a good averaging tool but close
        total_ts[i] = timeseries_results_inf
        total_ts_rec[i] = timeseries_results_rec

    infection_ts = np.mean(total_ts, axis=0)
    recover_ts = np.mean(total_ts_rec, axis=0)
    return ts, infection_ts, recover_ts
