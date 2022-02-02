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

LOCAL_DATA_PATH = "../data"

power_law_q2 = degree_distributions.power_law_degree_distrb(400, mu=10)
k_mean_degree = 2.5
# k_mean_degree = 5
er_degree_dist = degree_distributions.binomial_degree_distb(400, k_mean_degree)



# ensemble.run_ensemble_intervention_effects(er_degree_dist, '../data/rollouts/random/er_g5', 20000, 10000,
#                                            init_T=0.8, intervention_gen_list=[3, 4, 6], prop_reduced_list=[.2, .2, .3],
#                                            intervention_type="random_rollout", run_regular=False)

# ER with the magic random rollout effects
# pgf_formalism.compute_phase_space(20, 400, er_degree_dist, 0.8, True,
#                                   [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 16, 17, 18, 19],
#                                   f'{LOCAL_DATA_PATH}/rollouts/random/er_manatee_11_12_{0}',
#                                   rollout_dict={3: .2, 4: .2, 6: .3}, do_non_interv=False, do_interv=True,
#                                   intervention_type="random_rollout", pre_vax_correction=True)

plotting_util.plot_sims_vs_analytical_multigens([1, 2, 3, 4, 5, 6, 8], 200,
                                                f'{LOCAL_DATA_PATH}/rollouts/random/er_g5_intervene.txt', # non-intervention
                                                # Plotting intervention as the "no intervention" argument to just see those curves
                                                f'{LOCAL_DATA_PATH}/rollouts/random/er_manatee_11_12_{0}_intv.txt',
                                                same_plot=True, normalize_axis_x=False,
                                                colors_plot1={1: 'red', 2: 'orange', 3: 'blue', 4: 'purple', 5: 'pink',
                                                              6: 'grey', 8: 'darkgreen', 16: 'teal'})
# plt.savefig('random_rollout_ER_sample.png')
plt.title('Rollout at 3,4,6 of 20\%, 20\%, 30\%')
plt.savefig(f'{LOCAL_DATA_PATH}/manatee/rollouts/346_70pct_interventions.png')
plt.savefig(f'{LOCAL_DATA_PATH}/manatee/rollouts/346_70pct_interventions.svg', fmt='svg')
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