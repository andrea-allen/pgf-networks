import networkx as nx
import numpy as np
from epintervene.simobjects import simulation
from epintervene.simobjects import extended_simulation
from epintervene.simobjects import network
from src import pgf_formalism
import time
import matplotlib.pyplot as plt


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
    size_distrb_per_gen_intervention = simulate_ensemble(degree_distrb=degree_distrb, num_sims=num_sims, N=num_nodes,
                                                         intervention_gen=gen_intervene, intervention_T=T_intervene,
                                                         initial_T=init_T,
                                                         gamma=recover_rate, prop_reduced=prop_reduced,
                                                         intervention_gen_list=intervention_gen_list,
                                                         beta_redux_list=beta_redux_list,
                                                         prop_reduced_list=prop_reduced_list,
                                                         intervention_type=intervention_type
                                                         )
    print('Total time for ensemble was ', time.time() - start_time, ' with N=', num_nodes, ' nodes and ', num_sims,
          ' number of simulations.')
    np.savetxt(base_file_name + '_intervene.txt', size_distrb_per_gen_intervention,
               delimiter=',')


def run_ensemble_with_time_distribution(degree_distrb, base_file_name='sim_results',
                                        num_sims=10000, num_nodes=1000, init_T=0.8,
                                        gen_intervene=3, T_intervene=0.4, recover_rate=.001, prop_reduced=0.0,
                                        intervention_gen_list=None, beta_redux_list=None, prop_reduced_list=None,
                                        intervention_type="none",
                                        run_regular=True):
    start_time = time.time()
    t_buckets, time_buckets_distribution, generational_distribution, gen_emergence_avg = ensemble_time_distributions(degree_distrb, num_sims, num_nodes,
                                                                -1, 0.0, init_T, recover_rate, 0.0,
                                                                intervention_gen_list=None, beta_redux_list=None,
                                                                prop_reduced_list=None,
                                                                intervention_type="none")
    print('Total time for ensemble was ', time.time() - start_time, ' with N=', num_nodes, ' nodes and ', num_sims,
          ' number of simulations.')
    np.savetxt(base_file_name + '_generational.txt',
               generational_distribution,
               delimiter=',')
    np.savetxt(base_file_name + '_time_distribution.txt',
               time_buckets_distribution,
               delimiter=',')
    np.savetxt(base_file_name + '_time_distribution_values.txt',
               t_buckets,
               delimiter=',')
    np.savetxt(base_file_name + '_gen_emergence_times.txt',
               gen_emergence_avg,
               delimiter=',')


# Simulate ensemble saves results for generational time series, as opposed to real time results. Those are TBD
def simulate_ensemble(degree_distrb, num_sims=10, N=1000, intervention_gen=-1, intervention_T=0.0, initial_T=0.8,
                      gamma=0.1, prop_reduced=0.0, intervention_gen_list=None, beta_redux_list=None,
                      prop_reduced_list=None,
                      intervention_type="none"):
    # Configuring the parameters
    beta_init = -(gamma * initial_T) / (initial_T - 1)
    beta_interv = -(gamma * intervention_T) / (intervention_T - 1)

    # Setting up a results data structure
    outbreak_size_distrb_per_gen_matrix = np.zeros((100, N))

    # Generating the initial network
    G, pos = network.NetworkBuilder.from_degree_distribution(N, degree_distrb)
    A = np.array(nx.adjacency_matrix(G).todense())
    adjlist = network.NetworkBuilder.create_adjacency_list(G)
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
            adjlist = network.NetworkBuilder.create_adjacency_list(G)
            num_nodes_in_net = len(A[0])
            # Beta = np.full((num_nodes_in_net, num_nodes_in_net), beta_init)
            # Gamma = np.full(num_nodes_in_net, gamma)
            Beta = None
            Gamma = None
        # Get results of type generational time series
        results = run_single_simulation(A=A, adjlist=adjlist, Beta=Beta, Gamma=Gamma, current=i,
                                        results_type='generation', intervention_gen=intervention_gen,
                                        beta_interv=beta_interv, beta_init=beta_init, gamma_init=gamma,
                                        prop_reduced=prop_reduced,
                                        intervention_gen_list=intervention_gen_list, beta_redux_list=beta_redux_list,
                                        prop_reduced_list=prop_reduced_list, intervention_type=intervention_type)
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


def run_single_simulation(A, adjlist, Beta, Gamma, current, results_type='generation', intervention_gen=-1,
                          beta_interv=-1.0,
                          beta_init=1.0, gamma_init=0.001,
                          prop_reduced=0.0, intervention_gen_list=None, beta_redux_list=None, prop_reduced_list=None,
                          intervention_type="none"):
    start_time = time.time()
    N = len(A)
    # Constructing the simulation of specified Intervention Type, otherwise will run a regular simulation
    if intervention_type == "random-rollout":
        sim = extended_simulation.RandomRolloutSimulation(N=N, adjlist=adjlist, adjmatrix=A)
        sim.set_uniform_beta(beta_init)
        sim.set_uniform_gamma(gamma_init)
        sim.configure_intervention(intervention_gen_list=intervention_gen_list, beta_redux_list=beta_redux_list,
                                   proportion_reduced_list=prop_reduced_list)
        # TODO try this to see if can find bug:
        # sim.configure_intervention(intervention_gen_list=[5,6,7], beta_redux_list=[0,0,0],
        #                            proportion_reduced_list=[0,0,0])

    elif intervention_type == "random":
        sim = extended_simulation.RandomInterventionSim(N, adjmatrix=A, adjlist=adjlist)
        sim.set_uniform_beta(beta_init)
        sim.set_uniform_gamma(gamma_init)
        sim.configure_intervention(intervention_gen=intervention_gen, beta_redux=beta_interv,
                                   proportion_reduced=prop_reduced)
    elif intervention_type == "universal":
        sim = extended_simulation.UniversalInterventionSim(N, adjlist=adjlist)
        sim.set_uniform_beta(beta_init)
        sim.set_uniform_gamma(gamma_init)
        sim.configure_intervention(intervention_gen=intervention_gen, beta_redux=beta_interv)
    else:
        sim = simulation.Simulation(N=N, adj_list=adjlist, adj_matrix=A)
        sim.set_uniform_beta(beta_init)
        sim.set_uniform_gamma(gamma_init)

    # Run the simulation
    sim.run_sim(uniform_rate=True, wait_for_recovery=False)

    # Printing progress of the ensemble for the user to keep track
    if current % 5 == 0:
        print('current sim ' + str(current))
        print("--- %s seconds to run simulation---" % (time.time() - start_time))

    # Tabulating results based on user's specified input type
    if str(results_type).lower() == 'generation':
        generational_results = sim.tabulate_generation_results(100)
        return generational_results
    elif str(results_type).lower() == 'time':
        # TODO put the time buckets in here
        timeseries, timeseries_results_inf, timeseries_results_rec = sim.tabulate_continuous_time(1000,
                                                                                                  custom_range=True,
                                                                                                  custom_t_lim=20000)
        return timeseries, timeseries_results_inf, timeseries_results_rec
    elif str(results_type).lower() == 'time_and_gen':
        timeseries, timeseries_results_inf, timeseries_results_rec = sim.tabulate_continuous_time(1000,
                                                                                                  custom_range=True,
                                                                                                  custom_t_lim=20000)
        generational_results = sim.tabulate_generation_results(100)
        generational_emergence = sim.get_generational_emergence()
        return generational_results, generational_emergence, timeseries, timeseries_results_inf, timeseries_results_rec
    elif str(results_type).lower() == 'time_groups':
        timeseries, infection_timeseries_groups = sim.tabulate_continuous_time_with_groups(1000)
        return timeseries, infection_timeseries_groups
    else:
        print('No results type provided, returning basic time series results')
        timeseries, timeseries_results_inf, timeseries_results_rec = sim.tabulate_continuous_time(1000)

    return timeseries, timeseries_results_inf


#### In PROGRESS ####
def ensemble_time_distributions(degree_distrb, num_sims=10, N=1000, intervention_gen=-1, intervention_T=0.0,
                                initial_T=0.8,
                                gamma=0.1, prop_reduced=0.0, intervention_gen_list=None, beta_redux_list=None,
                                prop_reduced_list=None,
                                intervention_type="none"):
    # Configuring the parameters
    beta_init = -(gamma * initial_T) / (initial_T - 1)
    print(f'Beta {beta_init}')
    beta_interv = -(gamma * intervention_T) / (intervention_T - 1)

    # Setting up a results data structure
    ## This first definition of t_buckets does it in chunks of 1/beta
    t_buckets = np.arange(0, 15000, 1 / beta_init)
    ## This definition of t_buckets uses q*beta/(q*beta + gamma) *(q*beta + gamma).^{-1}, where q is avg excess degree
    q = pgf_formalism.z1_of(pgf_formalism.g1_of(degree_distrb))
    k = pgf_formalism.z1_of(degree_distrb)
    bucket_increment = (q*beta_init/(q*beta_init+gamma)) * (1/(q*beta_init + gamma))
    bucket_increment = (1/(q*beta_init + gamma)) # try T/(beta+gamma), then try again for a bigger q_bar
    bucket_increment = (beta_init/(beta_init+gamma)) / (beta_init+gamma)
    bucket_increment = 1/(q*(beta_init)) # Latest idea
    t_buckets = np.arange(0, 15000, bucket_increment)
    callibrated = False
    time_buckets_distribution = np.zeros((len(t_buckets), N))
    generational_distribution = np.zeros((100, N))
    # 2 rows: first is the average, 2nd is the number of times the generation actually ever existed (so the denominator of the average)
    gen_emergence_avg = np.zeros((2, 100))

    # Generating the initial network
    G, pos = network.NetworkBuilder.from_degree_distribution(N, degree_distrb)
    A = np.array(nx.adjacency_matrix(G).todense())
    adjlist = network.NetworkBuilder.create_adjacency_list(G)
    num_nodes_in_net = len(A[0])

    Beta = None
    Gamma = None


    # Running the ensemble
    for i in range(num_sims):
        # Generate a new network every 500 simulations
        if i % 500 == 0:
            G, pos = network.NetworkBuilder.from_degree_distribution(N, degree_distrb)
            A = np.array(nx.adjacency_matrix(G).todense())
            adjlist = network.NetworkBuilder.create_adjacency_list(G)
        # Get results of type generational time series and time series
        # generational_results, timeseries, timeseries_results_inf, timeseries_results_rec
        generational_results, generational_emergence, timeseries, timeseries_results_inf, timeseries_results_rec = \
            run_single_simulation(A=A, adjlist=adjlist, Beta=Beta, Gamma=Gamma, current=i, results_type='time_and_gen',
                                  intervention_gen=intervention_gen,
                                  beta_interv=beta_interv, beta_init=beta_init, gamma_init=gamma,
                                  prop_reduced=prop_reduced,
                                  intervention_gen_list=intervention_gen_list, beta_redux_list=beta_redux_list,
                                  prop_reduced_list=prop_reduced_list, intervention_type=intervention_type)
        # Recording ensemble results
        # Need to specify the range of gens based on beta
        if callibrated == False:
            new_time_buckets = np.zeros(len(t_buckets))
            time_values = timeseries
            curr_idx = 0
            curr_bucket_idx = 0
            for t in range(0, int(len(t_buckets))):
                found_nearest = False
                curr_bucket = t_buckets[curr_bucket_idx]
                while not found_nearest:
                    nearest = time_values[curr_idx]
                    if nearest >= curr_bucket:
                        new_time_buckets[curr_bucket_idx] = nearest
                        curr_bucket_idx += 1
                        found_nearest = True
                    else:
                        curr_idx += 1
            t_buckets = new_time_buckets
            callibrated = True

        # Recording results with time buckets
        # todo need the results to be cumulative
        # starting_tie = time.time()
        for t_time_idx in range(len(t_buckets)):
            # lazy way of cumulative infections is adding quantity recovered at that time too
            t_time = t_buckets[t_time_idx]
            t_time = int(t_time)
            t_idx = np.where(timeseries == t_time)
            num_total_infctd = int(timeseries_results_inf[t_idx] + timeseries_results_rec[t_idx])
            try:
                time_buckets_distribution[t_time_idx][num_total_infctd] += 1
            except IndexError:
                print('Index error for g: ', t_time_idx, ', gen_s: ', num_total_infctd)
                continue
        # print(f'Sorting the time took {time.time()-starting_tie}')

        # Recording ensemble results generational results
        for gen in range(len(generational_results)):
            num_total_infctd = int(generational_results[gen])
            try:
                generational_distribution[gen][num_total_infctd] += 1
            except IndexError:
                print('Index error for g: ', gen, ', gen_s: ', num_total_infctd)
                continue

        for gen in range(len(generational_results)):
            try:
                gen_emergence_avg[0][gen] += generational_emergence[gen]
                gen_emergence_avg[1][gen] += 1
            except KeyError:
                continue

        # Two other ideas:
        # 1. Average generational emergence, so a time series of the emergence of each generation (or a dict that's the average of each value)
        # 2. Time buckets, and average number of gens that have appeared by that time. Is this the same thing?
    # averaging results:
    for gen in range(100):
        gen_time_series = generational_distribution[gen]
        gen_time_series = gen_time_series / num_sims
        generational_distribution[gen] = gen_time_series
    # averaging results:
    for t in range(len(t_buckets)):
        t_time_series = time_buckets_distribution[t]
        t_time_series = t_time_series / num_sims
        time_buckets_distribution[t] = t_time_series
    # averaging results:
    for gen in range(100):
        if gen_emergence_avg[1][gen] > 0:
            gen_emergence_avg[0][gen] = gen_emergence_avg[0][gen] / gen_emergence_avg[1][gen]
    return t_buckets, time_buckets_distribution, generational_distribution, gen_emergence_avg[0]


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
        timeseries, timeseries_results_inf, timeseries_results_rec = run_single_simulation(A, Beta, Gamma, i,
                                                                                           'time_groups')
        ts = timeseries  # TODO is timeseries always the same? no, because it splits the max time into bins so maybe still not a good averaging tool but close
        total_ts[i] = timeseries_results_inf
        total_ts_rec[i] = timeseries_results_rec

    infection_ts = np.mean(total_ts, axis=0)
    recover_ts = np.mean(total_ts_rec, axis=0)
    return ts, infection_ts, recover_ts
