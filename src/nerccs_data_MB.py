
from src import degree_distributions

from src import pgf_formalism
from analysis import ensemble



LOCAL_DATA_PATH = "../data"
max_len_rand = 7000
max_len_targ = 7000
gens_compute = 13
power_law_q2 = degree_distributions.power_law_degree_distrb(400, mu=10)
k_mean_degree = 2.5
# k_mean_degree = 5
er_degree_dist = degree_distributions.binomial_degree_distb(400, k_mean_degree)

geom_degree_dist_rand = degree_distributions.exponential_degree_dist(max_len_rand, .6) # n and then p
geom_degree_dist_targ = degree_distributions.exponential_degree_dist(max_len_targ, .6) # n and then p
"""
EXAMPLES
###########
#targeted intervention simulations on powerlaw
## SET UP YOUR OWN FILE PATHS BELOW
"""

# Non intervention
# pgf_formalism.compute_phase_space(gens_compute, max_len_targ, geom_degree_dist_targ,  0.8, True,
#                                   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12],
#                                   f'../nerccs2022/geo_05-5-22_non_interv_{{0}}', do_non_interv=True, do_interv=False, pre_vax_correction=True)
#

#Random (single and roll outs)

pgf_formalism.compute_phase_space(gens_compute, max_len_rand, geom_degree_dist_rand,  0.8, True,
                                  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12],
                                  f'../nerccs2022/geo_05-11-22_random_gen5_15_{{0}}',
                                  rollout_dict={5:.15}, do_non_interv=False, do_interv=True,
                                  intervention_type="random_rollout", pre_vax_correction=True)

# ensemble.run_ensemble_intervention_effects(geom_degree_dist_rand, f'../nerccs2022/geo_05-5-22_random_gen5_sims', num_sims=75000, num_nodes=20000,
#                                            init_T=0.8, intervention_gen_list=[5], prop_reduced_list=[.15], intervention_type="random_rollout", run_regular=False)
#
#
pgf_formalism.compute_phase_space(gens_compute, max_len_rand, geom_degree_dist_rand, 0.8, True,
                                    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12],
                                  f'../nerccs2022/geo_05-11-22_random_rollout_10_{{0}}',
                                  rollout_dict={4:0.1, 6:0.1, 8:0.1, 10:0.1, 12:0.1}, do_non_interv=False, do_interv=True,
                                  intervention_type="random_rollout", pre_vax_correction=True)

# ensemble.run_ensemble_intervention_effects(geom_degree_dist_rand, f'../nerccs2022/geo_05-5-22_random_rollout_15_sims', num_sims=75000, num_nodes=20000,
#                                            init_T=0.8, intervention_gen_list=[4, 6, 8, 10, 12], prop_reduced_list=[0.1, 0.1, 0.1, 0.1, 0.1], intervention_type="random_rollout", run_regular=False)
#



# Targeted (single and roll outs)
# Do this
pgf_formalism.compute_phase_space(gens_compute, max_len_targ, geom_degree_dist_targ,  0.8, True,
                                  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                                  f'../nerccs2022/geo_05-11-22_targeted_gen5_10_{{0}}',
                                  rollout_dict={5:.1}, do_non_interv=False, do_interv=True,
                                  intervention_type="targeted_rollout", pre_vax_correction=True)

# ensemble.run_ensemble_intervention_effects(geom_degree_dist_targ, f'../nerccs2022/geo_05-10-22_targeted_gen5_sims_2p5', num_sims=75000, num_nodes=20000,
#                                            init_T=0.8, intervention_gen_list=[5], prop_reduced_list=[.025], intervention_type="targeted_rollout", run_regular=False)
#
#Do this
pgf_formalism.compute_phase_space(gens_compute, max_len_targ, geom_degree_dist_targ,  0.8, True,
                                  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                                  f'../nerccs2022/geo_05-11-22_targeted_rollout_2_{{0}}',
                                  rollout_dict={4:.02, 6:.02, 8:.02, 10:.02, 12:.02}, do_non_interv=False, do_interv=True,
                                  intervention_type="targeted_rollout", pre_vax_correction=True)

# ensemble.run_ensemble_intervention_effects(geom_degree_dist_targ, f'../nerccs2022/geo_05-4-22_targeted_rollout_2_sims', num_sims=75000, num_nodes=20000,
#                                            init_T=0.8, intervention_gen_list=[4,6,8,10,12], prop_reduced_list=[.02,.02,.02,.02,.02],
#                                          intervention_type="targeted_rollout", run_regular=False)

# DO This
# pgf_formalism.compute_phase_space(13, max, geom_degree_dist,  0.8, True,
#                                   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19],
#                                   f'../nerccs2022/geo_05-11-22_Non_interv_{{0}}', do_non_interv=True, do_interv=False,
#                                   pre_vax_correction=True)



# pgf_formalism.compute_phase_space(20, 400, geom_degree_dist,  0.8, True,
#                                   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19],
#                                   f'../nerccs2022/geo_04-25-22_targeted_rollout_2point5per_{{0}}',
#                                   rollout_dict={4:.025, 6:.025, 8:.025, 10:.025, 12:.025}, do_non_interv=False, do_interv=True,
#                                   intervention_type="targeted_rollout", pre_vax_correction=True)
#
# pgf_formalism.compute_phase_space(20, 400, geom_degree_dist,  0.8, True,
#                                   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12],
#                                   f'../nerccs2022/geo_04-25-22_targeted_rollout_5per_{{0}}',
#                                   rollout_dict={4:.05, 6:.05, 8:.05, 10:.05, 12:.05}, do_non_interv=False, do_interv=True,
#                                   intervention_type="targeted_rollout", pre_vax_correction=True)
#
#
# pgf_formalism.compute_phase_space(20, 400, geom_degree_dist,  0.8, True,
#                                   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12],
#                                   f'../nerccs2022/geo_04-25-22_targeted_rollout_10per_{{0}}',
#                                   rollout_dict={4:.1, 6:.1, 8:.1, 10:.1, 12:.1}, do_non_interv=False, do_interv=True,
#                                   intervention_type="targeted_rollout", pre_vax_correction=True)

# pgf_formalism.compute_phase_space(20, 400, geom_degree_dist,  0.8, True,
#                                   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12],
#                                   f'../nerccs2022/geo_04-25-22_targeted_rollout_20per_{{0}}',
#                                   rollout_dict={4:.2, 6:.2, 8:.2, 10:.2, 12:.2}, do_non_interv=False, do_interv=True,
#                                   intervention_type="targeted_rollout", pre_vax_correction=True)

#
# for decimal, percent in zip([0.01, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.20],[2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20]):
#     pgf_formalism.compute_phase_space(gens_compute, max_len_targ, geom_degree_dist_targ,  0.8, True,
#                                       [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12],
#                                       f'../nerccs2022/geo_05-5-22_targeted_rollout_{percent}_resultsfigs_{{0}}',
#                                       rollout_dict={4:decimal, 6:decimal, 8:decimal, 10:decimal, 12:decimal}, do_non_interv=False, do_interv=True,
#                                       intervention_type="targeted_rollout", pre_vax_correction=True)

# pgf_formalism.compute_phase_space(gens_compute, max_len_targ, geom_degree_dist_targ, 0.8, True,
#                                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12],
#                                     f'../nerccs2022/geo_05-5-22_targeted_rollout_1_resultsfigs_{{0}}',
#                                     rollout_dict={4: 0.01, 6: 0.01, 8: 0.01, 10:0.01, 12: 0.01},
#                                     do_non_interv=False, do_interv=True,
#                                     intervention_type="targeted_rollout", pre_vax_correction=True)