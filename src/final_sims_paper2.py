from src import degree_distributions

from src import pgf_formalism
from analysis import ensemble



LOCAL_DATA_PATH = "../data"
max_len_rand = 7000
max_len_targ = 7000
gens_compute = 16
power_law_q2 = degree_distributions.power_law_degree_distrb(400, mu=10)
k_mean_degree = 2.5
# k_mean_degree = 5
er_degree_dist = degree_distributions.binomial_degree_distb(400, k_mean_degree)

geom_degree_dist_rand = degree_distributions.exponential_degree_dist(max_len_rand, .6) # n and then p
geom_degree_dist_targ = degree_distributions.exponential_degree_dist(max_len_targ, .6) # n and then p

# Random at 0.15
pgf_formalism.compute_phase_space(gens_compute, max_len_rand, geom_degree_dist_rand,  0.8, True,
                                  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16],
                                  f'../nerccs2022/geo_10-26-22_random_gen5_point15_{{0}}',
                                  rollout_dict={5:.015}, do_non_interv=False, do_interv=True,
                                  intervention_type="random_rollout", pre_vax_correction=True)

ensemble.run_ensemble_intervention_effects(geom_degree_dist_rand, f'../nerccs2022/geo_10-26-22_random_gen5_sims', num_sims=75000, num_nodes=20000,
                                           init_T=0.8, intervention_gen_list=[5], prop_reduced_list=[.015], intervention_type="random_rollout", run_regular=False)

#Random rollout with 0.25 cumulative intervention

pgf_formalism.compute_phase_space(gens_compute, max_len_rand, geom_degree_dist_rand, 0.8, True,
                                    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16],
                                  f'../nerccs2022/geo_10-26-22_random_rollout_005_{{0}}',
                                  rollout_dict={4:0.005, 6:0.005, 8:0.005, 10:0.005, 12:0.005}, do_non_interv=False, do_interv=True,
                                  intervention_type="random_rollout", pre_vax_correction=True)

ensemble.run_ensemble_intervention_effects(geom_degree_dist_rand, f'../nerccs2022/geo_10-26-22_random_rollout_005_sims', num_sims=75000, num_nodes=20000,
                                           init_T=0.8, intervention_gen_list=[4, 6, 8, 10, 12], prop_reduced_list=[0.005, 0.005, 0.005, 0.005, 0.005], intervention_type="random_rollout", run_regular=False)


#comparison figure
for decimal, percent in zip([0.001, 0.0025, 0.005, 0.0075, 0.010, 0.0125, 0.015, 0.0175, 0.020],[0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.20]):
    pgf_formalism.compute_phase_space(gens_compute, max_len_targ, geom_degree_dist_targ,  0.8, True,
                                      [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16],
                                      f'../nerccs2022/geo_10-26-22_targeted_rollout_{percent}_resultsfigs_{{0}}',
                                      rollout_dict={4:decimal, 6:decimal, 8:decimal, 10:decimal, 12:decimal}, do_non_interv=False, do_interv=True,
                                      intervention_type="targeted_rollout", pre_vax_correction=True)
