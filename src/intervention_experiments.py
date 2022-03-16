"""
File to try out various intervention modeling methods, run and compare against simulations, general scratchwork
"""
## Goal today (1/25): Not to figure out what's wrong with manatee, just make the code more manageable.

import matplotlib.pyplot as plt
from matplotlib import rc
from src import degree_distributions
from src import plotting_util
from src import figs_for_paper
from analysis import covid
import numpy as np
import random
from src import pgf_formalism
from analysis import ensemble
import src.manatee as mnte
import time

LOCAL_DATA_PATH = "../data"

power_law_q2 = degree_distributions.power_law_degree_distrb(400, mu=10)
k_mean_degree = 2.5
# k_mean_degree = 5
er_degree_dist = degree_distributions.binomial_degree_distb(400, k_mean_degree)

# plt.figure('4-6-8 90% POWER LAW')
# plotting_util.plot_sims_vs_analytical_multigens([1, 2, 3, 4, 5, 9, 12, 13, 15, 18], 200,
#                                                 f'{LOCAL_DATA_PATH}/manatee/rollouts/pl_02-22-22_g4-6-8_intervene.txt', # non-intervention
#                                                 # Plotting intervention as the "no intervention" argument to just see those curves
#                                                 # f'{LOCAL_DATA_PATH}/rollouts/random/er_manatee_11_12_{{0}}_intv.txt',
#                                                 f'{LOCAL_DATA_PATH}/manatee/rollouts/pl_02-22-22_g4-6-8_{{0}}_intv.txt',
#                                                 same_plot=True, normalize_axis_x=False,
#                                                 colors_plot1={1: 'red', 2: 'orange', 3: 'brown', 4: 'green', 5: 'blue',
#                                                               6: 'purple', 9: 'pink', 12: 'teal', 13: 'black', 15:'magenta', 18:'grey'})
# plt.title('4-6-8 90% POWER LAW')
# plt.show()

# plotting_util.plot_sims_vs_analytical_multigens([1, 2, 3, 4, 5, 6, 8], 200,
#                                                 f'{LOCAL_DATA_PATH}/rollouts/random/er_g5_intervene.txt', # non-intervention
#                                                 # Plotting intervention as the "no intervention" argument to just see those curves
#                                                 # f'{LOCAL_DATA_PATH}/rollouts/random/er_manatee_11_12_{{0}}_intv.txt',
#                                                 f'{LOCAL_DATA_PATH}/rollouts/random/er_manatee_02_15_alt3_{{0}}_intv.txt',
#                                                 same_plot=True, normalize_axis_x=False,
#                                                 colors_plot1={1: 'red', 2: 'orange', 3: 'blue', 4: 'purple', 5: 'pink',
#                                                               6: 'grey', 8: 'darkgreen', 15: 'teal', 18: 'black'})
# # plt.savefig('random_rollout_ER_sample.png')
# plt.title('Rollout at 3,4,6 of 20\%, 20\%, 30\%')
# # plt.savefig(f'{LOCAL_DATA_PATH}/manatee/rollouts/346_70pct_interventions.png')
# # plt.savefig(f'{LOCAL_DATA_PATH}/manatee/rollouts/346_70pct_interventions.svg', fmt='svg')
# plt.show()

############
# targeted intervention simulations on powerlaw
# start_time = time.time()
# ensemble.run_ensemble_intervention_effects(power_law_q2, '../data/manatee/rollouts/pl_03-16-22_targeted_6pct', 10000, 10000,
#                                            init_T=0.8, intervention_gen_list=[4,6,8], prop_reduced_list=[.02, .02, .02],
#                                            intervention_type="targeted_rollout", run_regular=False)
# # end_time = time.time() - start_time
# # np.savetxt(f'{LOCAL_DATA_PATH}/manatee/rollouts/plaw_3-16_targeted_time.txt', np.array([end_time]))
# pgf_formalism.compute_phase_space(20, 400, power_law_q2,  0.8, True,
#                                   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 17, 18, 19],
#                                   f'{LOCAL_DATA_PATH}/manatee/rollouts/pl_03-16-22_targeted_6pct_{{0}}',
#                                   rollout_dict={4:.02, 6:.02, 8:.02}, do_non_interv=False, do_interv=True,
#                                   intervention_type="targeted_rollout", pre_vax_correction=True)
#
# plotting_util.plot_sims_vs_analytical_multigens([1, 2, 3, 4, 5, 9, 12, 13, 15, 18], 200,
#                                                 f'{LOCAL_DATA_PATH}/manatee/rollouts/pl_03-16-22_targeted_6pct_intervene.txt', # non-intervention
#                                                 # Plotting intervention as the "no intervention" argument to just see those curves
#                                                 # f'{LOCAL_DATA_PATH}/rollouts/random/er_manatee_11_12_{{0}}_intv.txt',
#                                                 f'{LOCAL_DATA_PATH}/manatee/rollouts/pl_03-16-22_targeted_6pct_{{0}}_intv.txt',
#                                                 same_plot=True, normalize_axis_x=False,
#                                                 colors_plot1={1: 'red', 2: 'orange', 3: 'brown', 4: 'green', 5: 'blue',
#                                                               6: 'purple', 9: 'pink', 12: 'teal', 13: 'black', 15:'magenta', 18:'grey'})
# # plt.title('SIMS ONLY power law targeted')
# plt.show()

###########
# targeted intervention simulations on powerlaw
# start_time = time.time()
# CODE TO SIMULATE
# ensemble.run_ensemble_intervention_effects(power_law_q2, '../data/manatee/rollouts/pl_03-16-22_targeted_132pct', 10000, 10000,
#                                            init_T=0.8, intervention_gen_list=[4,6,8], prop_reduced_list=[.01, .03, .02],
#                                            intervention_type="targeted_rollout", run_regular=False)
# end_time = time.time() - start_time
# np.savetxt(f'{LOCAL_DATA_PATH}/manatee/rollouts/plaw_3-16_targeted_time.txt', np.array([end_time]))
# CODE TO MODEL
# pgf_formalism.compute_phase_space(20, 400, power_law_q2,  0.8, True,
#                                   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 17, 18, 19],
#                                   f'{LOCAL_DATA_PATH}/manatee/rollouts/pl_03-16-22_targeted_132pct_{{0}}',
#                                   rollout_dict={4:.01, 6:.03, 8:.02}, do_non_interv=False, do_interv=True,
#                                   intervention_type="targeted_rollout", pre_vax_correction=True)

# CODE TO PLOT
plotting_util.plot_sims_vs_analytical_multigens([1, 2, 3, 4, 5, 9, 12, 13, 15, 18], 200,
                                                f'{LOCAL_DATA_PATH}/manatee/rollouts/pl_03-16-22_targeted_132pct_intervene.txt', # non-intervention
                                                # Plotting intervention as the "no intervention" argument to just see those curves
                                                # f'{LOCAL_DATA_PATH}/rollouts/random/er_manatee_11_12_{{0}}_intv.txt',
                                                f'{LOCAL_DATA_PATH}/manatee/rollouts/pl_03-16-22_targeted_132pct_{{0}}_intv.txt',
                                                same_plot=True, normalize_axis_x=False,
                                                colors_plot1={1: 'red', 2: 'orange', 3: 'brown', 4: 'green', 5: 'blue',
                                                              6: 'purple', 9: 'pink', 12: 'teal', 13: 'black', 15:'magenta', 18:'grey'})
# plt.title('SIMS ONLY power law targeted')
plt.show()



############


############
# targeted intervention simulations on erdos-renyi
start_time = time.time()
# ensemble.run_ensemble_intervention_effects(er_degree_dist, '../data/manatee/rollouts/er_03-16-22_targeted_12pct_3rolls', 3000, 10000,
#                                            init_T=0.8, intervention_gen_list=[4,6, 8], prop_reduced_list=[.04, .04, .04],
#                                            intervention_type="targeted_rollout", run_regular=False)
# end_time = time.time() - start_time
# np.savetxt(f'{LOCAL_DATA_PATH}/manatee/rollouts/er_3-08_targeted_time.txt', np.array([end_time]))
# pgf_formalism.compute_phase_space(20, 400, er_degree_dist,  0.8, True,
#                                   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 17, 18, 19],
#                                   f'{LOCAL_DATA_PATH}/manatee/rollouts/er_03-16-22_targeted_12pct_3rolls_{{0}}',
#                                   rollout_dict={4:.04, 6:.04, 8:.04}, do_non_interv=False, do_interv=True,
#                                   intervention_type="targeted_rollout", pre_vax_correction=True)

plotting_util.plot_sims_vs_analytical_multigens([1, 2, 3, 4, 5, 9, 12, 13, 15, 18], 200,
                                                f'{LOCAL_DATA_PATH}/manatee/rollouts/er_03-16-22_targeted_12pct_3rolls_intervene.txt', # non-intervention
                                                # Plotting intervention as the "no intervention" argument to just see those curves
                                                # f'{LOCAL_DATA_PATH}/rollouts/random/er_manatee_11_12_{{0}}_intv.txt',
                                                f'{LOCAL_DATA_PATH}/manatee/rollouts/er_03-16-22_targeted_12pct_3rolls_{{0}}_intv.txt',
                                                same_plot=True, normalize_axis_x=False,
                                                colors_plot1={1: 'red', 2: 'orange', 3: 'brown', 4: 'green', 5: 'blue',
                                                              6: 'purple', 9: 'pink', 12: 'teal', 13: 'black', 15:'magenta', 18:'grey'})
# plt.title('ER targeted 6% 2 rolls')
plt.title('ER targeted 12% 3 rolls')
plt.show()

############
# new examples trying to completely kill off the epidemic for proof of concept
start_time = time.time()
# ensemble.run_ensemble_intervention_effects(power_law_q2, '../data/manatee/rollouts/pl_02-22-22_g4-6-8', 2000, 10000,
#                                            init_T=0.8, intervention_gen_list=[4,6,8], prop_reduced_list=[.3, .3, .3],
#                                            intervention_type="random_rollout", run_regular=False)
end_time = time.time() - start_time
# np.savetxt(f'{LOCAL_DATA_PATH}/manatee/rollouts/plaw_2-22_time.txt', np.array([end_time]))

# pgf_formalism.compute_phase_space(20, 400, power_law_q2,  0.8, True,
#                                   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 17, 18, 19],
#                                   f'{LOCAL_DATA_PATH}/manatee/rollouts/pl_02-22-22_g4-6-8_{{0}}',
#                                   rollout_dict={4:.3, 6:.3, 8:.3}, do_non_interv=False, do_interv=True,
#                                   intervention_type="random_rollout", pre_vax_correction=True)
plt.figure('4-6-8 90% POWER LAW')
plotting_util.plot_sims_vs_analytical_multigens([1, 2, 3, 4, 5, 9, 12, 13, 15, 18], 200,
                                                f'{LOCAL_DATA_PATH}/manatee/rollouts/pl_02-22-22_g4-6-8_intervene.txt', # non-intervention
                                                # Plotting intervention as the "no intervention" argument to just see those curves
                                                # f'{LOCAL_DATA_PATH}/rollouts/random/er_manatee_11_12_{{0}}_intv.txt',
                                                f'{LOCAL_DATA_PATH}/manatee/rollouts/pl_02-22-22_g4-6-8_{{0}}_intv.txt',
                                                same_plot=True, normalize_axis_x=False,
                                                colors_plot1={1: 'red', 2: 'orange', 3: 'brown', 4: 'green', 5: 'blue',
                                                              6: 'purple', 9: 'pink', 12: 'teal', 13: 'black', 15:'magenta', 18:'grey'})
plt.title('4-6-8 90% POWER LAW')
# plt.show()

############
# new examples with a slower rollout for random vaccination, with powerlaw
# start_time = time.time()
# ensemble.run_ensemble_intervention_effects(power_law_q2, '../data/manatee/rollouts/pl_02-22-22_g5-6-7-8-9', 10000, 10000,
#                                            init_T=0.8, intervention_gen_list=[5, 6, 7, 8, 9], prop_reduced_list=[.1, .1, .1, .1, .1],
#                                            intervention_type="random_rollout", run_regular=False)
# end_time = time.time() - start_time
# np.savetxt(f'{LOCAL_DATA_PATH}/manatee/rollouts/plaw_2-22_time.txt', np.array([end_time]))
#
# pgf_formalism.compute_phase_space(20, 400, power_law_q2,  0.8, True,
#                                   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 17, 18, 19],
#                                   f'{LOCAL_DATA_PATH}/manatee/rollouts/pl_02-22-22_g5-6-7-8-9_{{0}}',
#                                   rollout_dict={5:.1, 6:.1, 7:.1, 8:.1, 9:.1}, do_non_interv=False, do_interv=True,
#                                   intervention_type="random_rollout", pre_vax_correction=True)
plt.figure('5-6-7-8-9 POWER LAW')
plotting_util.plot_sims_vs_analytical_multigens([1, 2, 3, 4, 5, 9, 12, 13, 15, 18], 200,
                                                f'{LOCAL_DATA_PATH}/manatee/rollouts/pl_02-22-22_g5-6-7-8-9_intervene.txt', # non-intervention
                                                # Plotting intervention as the "no intervention" argument to just see those curves
                                                # f'{LOCAL_DATA_PATH}/rollouts/random/er_manatee_11_12_{{0}}_intv.txt',
                                                f'{LOCAL_DATA_PATH}/manatee/rollouts/pl_02-22-22_g5-6-7-8-9_{{0}}_intv.txt',
                                                same_plot=True, normalize_axis_x=False,
                                                colors_plot1={1: 'red', 2: 'orange', 3: 'brown', 4: 'green', 5: 'blue',
                                                              6: 'purple', 9: 'pink', 12: 'teal', 13: 'black', 15:'magenta', 18:'grey'})
plt.title('5-6-7-8-9 POWER LAW')
# plt.show()

##################
# start_time = time.time()
# ensemble.run_ensemble_intervention_effects(er_degree_dist, '../data/manatee/rollouts/er_02-22-22_g5-6-7-8-9', 10000, 10000,
#                                            init_T=0.8, intervention_gen_list=[5, 6, 7, 8, 9], prop_reduced_list=[.1, .1, .1, .1, .1],
#                                            intervention_type="random_rollout", run_regular=False)
# end_time = time.time() - start_time
# np.savetxt(f'{LOCAL_DATA_PATH}/manatee/rollouts/ER_2-22_time.txt', np.array([end_time]))
#
# pgf_formalism.compute_phase_space(20, 400, er_degree_dist,  0.8, True,
#                                   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 17, 18, 19],
#                                   f'{LOCAL_DATA_PATH}/manatee/rollouts/er_02-22-22_g5-6-7-8-9_{{0}}',
#                                   rollout_dict={5:.1, 6:.1, 7:.1, 8:.1, 9:.1}, do_non_interv=False, do_interv=True,
#                                   intervention_type="random_rollout", pre_vax_correction=True)
plt.figure('5-6-7-8-9 ERDOS RENYI')
plotting_util.plot_sims_vs_analytical_multigens([1, 2, 3, 4, 5, 9, 12, 13, 15, 18], 200,
                                                f'{LOCAL_DATA_PATH}/manatee/rollouts/er_02-22-22_g5-6-7-8-9_intervene.txt', # non-intervention
                                                # Plotting intervention as the "no intervention" argument to just see those curves
                                                # f'{LOCAL_DATA_PATH}/rollouts/random/er_manatee_11_12_{{0}}_intv.txt',
                                                f'{LOCAL_DATA_PATH}/manatee/rollouts/er_02-22-22_g5-6-7-8-9_{{0}}_intv.txt',
                                                same_plot=True, normalize_axis_x=False,
                                                colors_plot1={1: 'red', 2: 'orange', 3: 'brown', 4: 'green', 5: 'blue',
                                                              6: 'purple', 9: 'pink', 12: 'teal', 13: 'black', 15:'magenta', 18:'grey'})
plt.title('5-6-7-8-9 ERDOS RENYI')
# plt.show()
###################




# ensemble.run_ensemble_intervention_effects(power_law_q2, '../data/rollouts/random/pl_02-22-22_g5-9-13', 6000, 10000,
#                                            init_T=0.8, intervention_gen_list=[5, 9, 13], prop_reduced_list=[.2, .2, .1],
#                                            intervention_type="random_rollout", run_regular=False)
#
# pgf_formalism.compute_phase_space(20, 400, power_law_q2,  0.8, True,
#                                   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 17, 18, 19],
#                                   f'{LOCAL_DATA_PATH}/rollouts/random/pl_02-22-22_g5-9-13_{{0}}',
#                                   rollout_dict={5: .2, 9: .2, 13: .1}, do_non_interv=False, do_interv=True,
#                                   intervention_type="random_rollout", pre_vax_correction=True)

plt.figure('5-9-13 POWER LAW')
plotting_util.plot_sims_vs_analytical_multigens([1, 2, 3, 4, 5, 9, 12, 13, 15, 18], 200,
                                                f'{LOCAL_DATA_PATH}/rollouts/random/pl_02-22-22_g5-9-13_intervene.txt', # non-intervention
                                                # Plotting intervention as the "no intervention" argument to just see those curves
                                                # f'{LOCAL_DATA_PATH}/rollouts/random/er_manatee_11_12_{{0}}_intv.txt',
                                                f'{LOCAL_DATA_PATH}/rollouts/random/pl_02-22-22_g5-9-13_{{0}}_intv.txt',
                                                same_plot=True, normalize_axis_x=False,
                                                colors_plot1={1: 'red', 2: 'orange', 3: 'brown', 4: 'green', 5: 'blue',
                                                              6: 'purple', 9: 'pink', 12: 'teal', 13: 'black', 15:'magenta', 18:'grey'})
# plt.savefig('random_rollout_ER_sample.png')
plt.title('Slow rollout at 5,9,13 of 20\%, 20\%, 10\%')
# plt.savefig(f'{LOCAL_DATA_PATH}/manatee/rollouts/346_70pct_interventions.png')
# plt.savefig(f'{LOCAL_DATA_PATH}/manatee/rollouts/346_70pct_interventions.svg', fmt='svg')
# plt.show()

############



# ensemble.run_ensemble_intervention_effects(er_degree_dist, '../data/rollouts/random/er_g5', 20000, 10000,
#                                            init_T=0.8, intervention_gen_list=[3, 4, 6], prop_reduced_list=[.2, .2, .3],
#                                            intervention_type="random_rollout", run_regular=False)

# ER with the magic random rollout effects
# pgf_formalism.compute_phase_space(20, 400, er_degree_dist, 0.8, True,
#                                   [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 16, 17, 18, 19],
#                                   f'{LOCAL_DATA_PATH}/rollouts/random/er_manatee_02_15_alt3_{{0}}',
#                                   rollout_dict={3: .2, 4: .2, 6: .3}, do_non_interv=False, do_interv=True,
#                                   intervention_type="random_rollout", pre_vax_correction=True)

plotting_util.plot_sims_vs_analytical_multigens([1, 2, 3, 4, 5, 6, 8], 200,
                                                f'{LOCAL_DATA_PATH}/rollouts/random/er_g5_intervene.txt', # non-intervention
                                                # Plotting intervention as the "no intervention" argument to just see those curves
                                                # f'{LOCAL_DATA_PATH}/rollouts/random/er_manatee_11_12_{{0}}_intv.txt',
                                                f'{LOCAL_DATA_PATH}/rollouts/random/er_manatee_02_15_alt3_{{0}}_intv.txt',
                                                same_plot=True, normalize_axis_x=False,
                                                colors_plot1={1: 'red', 2: 'orange', 3: 'blue', 4: 'purple', 5: 'pink',
                                                              6: 'grey', 8: 'darkgreen', 15: 'teal', 18: 'black'})
# plt.savefig('random_rollout_ER_sample.png')
plt.title('Rollout at 3,4,6 of 20\%, 20\%, 30\%')
# plt.savefig(f'{LOCAL_DATA_PATH}/manatee/rollouts/346_70pct_interventions.png')
# plt.savefig(f'{LOCAL_DATA_PATH}/manatee/rollouts/346_70pct_interventions.svg', fmt='svg')
plt.show()

""" An example of a single-batch intervention that works with current manatee (data is local)"""
plotting_util.plot_sims_vs_analytical_multigens([1, 3, 5, 6, 8, 10], 400,
                                                # None,
                                                # '../data/paper_2_figs_results/random/er_g5_pois2{0}_intv.txt',
                                                f'{LOCAL_DATA_PATH}/manatee/random/er_g5_40pct_intervene.txt', #Plotting intervention as the "no intervention" argument to just see those curves
                                                f'{LOCAL_DATA_PATH}/manatee/random/er_g5_40pct_{0}_intv.txt',
                                                same_plot=True, normalize_axis_x=False,
                                                colors_plot1={1: 'darkblue', 3: 'mediumblue', 5: 'blue',
                                                              6: 'dodgerblue', 8: 'cadetblue', 10:'cornflowerblue'},
                                                colors_plot2={1: 'indigo', 3: 'darkorchid', 5: 'mediumpurple',
                                                              6: 'hotpink', 8: 'plum', 10:'darkblue'})
plt.title('Intervention at Generation 5 of 40% random vaccination')
# plt.savefig(f'{LOCAL_DATA_PATH}/manatee_single_interv_figure.png')
# plt.savefig(f'{LOCAL_DATA_PATH}/manatee_single_interv_figure.svg', fmt='svg')
plt.show()