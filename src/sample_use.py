import time

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy import stats

from analysis import ensemble
from src import degree_distributions
from src import pgf_formalism
from src import plotting_util
from src import figs_for_paper
from analysis import covid

import epintervene

# ## CREATE A DEGREE DISTRIBUTION
power_law_dd = degree_distributions.power_law_degree_distrb(400)
power_law_q3 = degree_distributions.power_law_degree_distrb(2000, mu=10)  # q is 3, k is 1.7
poisson = degree_distributions.binomial_degree_distb(2000, 2.5)

# ## COMPUTE AND SAVE THE PHASE SPACE MATRICES FOR SPECIFIED GENERATIONS
pgf_formalism.compute_phase_space(num_gens=17, num_nodes=2000, degree_distribution=poisson, transmissibility=0.8,
                                  save_results=True, gens_to_save=[2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 16],
                                  file_fmt_to_save='../data/paper/ER_q2.5_T8_ifft_g{0}',
                                  do_non_interv=True, do_interv=False)

# ## PLOT SIMULATION RESULTS AGAINST THEORETICAL RESULTS, BOTH READ FROM FILES: (LOOK INSIDE FUNCTION FOR FILE FORMAT
# AND LOCATION SPECIFICATION)
plotting_util.plot_sims_vs_analytical_multigens([3, 6, 10, 12], 400,
                                                f'../data/testing/plaw_T8_10k_80ksims_generational.txt',
                                                '../data/testing/powerlaw_T8_g{0}.txt',
                                                same_plot=True, normalize_axis_x=False)
plt.show()

# ## RUN SPECIFIC EXPERIMENTS FOR OUR PAPER
# RUNS BOTH SIMULATIONS AND THEORETICAL COMPUTATIONS, SAVES FILES
# CHECK CODE FOR FILE FMT SPECIFICATION AND LOCATION
figs_for_paper.erdos_renyi_exp(file_root='erdos_renyi_10k_10ksims_pt1', num_sims=10000, num_nodes=10000, kill_by=16)
# ACTIVE GEN SIZES WILL RECORD ACTIVE GENS DURING SIMULATION
# THIS SIGNIFICANTLY SLOWS DOWN SIMS BUT REQUIRES FEWER TO GENERATE SMOOTH RESULTS
figs_for_paper.power_law_exp(file_root='plaw_T8_10k_1ksims_q3_p7', num_sims=1000, num_nodes=10000, active_gen_sizes_on=True)
# PLOTTING FUNCTION FOR FIGURES FOR OUR PAPER, PLUS ADDITIONAL UNUSED FIGURES/ANALYSIS
figs_for_paper.results_plots(file_root='plaw_T8_10k_40ksims_q3_p6', q_degree=3.04, active_gen_sizes_on=False)

# ## PLOTTING PHASE SPACE FROM DATA
# plotting_util.phaseSpace_from_data('../data/talk_results/powerlawrollout4565.txt', 5, '$(s,m)$ phase space')

# ## RUNNING ENSEMBLE FOR INTERVENTION EFFECTS
### Universal intervention:
ensemble.run_ensemble_intervention_effects(power_law_dd, '../data/talk_results/pl_T8', 75000, 10000,
                                           init_T=0.8, gen_intervene=4, T_intervene=0.4,
                                        intervention_type="none", run_regular=True)

# ## PLOTS FOR NERCCS TALK, USED FOR INTERVENTION THEORY
# plotting_util.plots_for_nerccs_talk([3, 6, 10, 15], 400,  '../data/talk_results/pl_T8_no_intervene.txt', '../data/talk_results/powerlawrollout456{0}.txt',
#                                                 '../data/talk_results/universal_pl_T8_gen4_T4_intervene.txt', '../data/talk_results/powerlawuniversal4_T4{0}_intv.txt',
#                                                 same_plot=True, normalize_axis_x=False, plot_distribution_inset=True, inset_label=f'$T_g \\rightarrow T_4=0.4$')
#
# ## Universal Intervention
# # plotting_util.plots_for_nerccs_talk([3, 6, 10, 15], 400,  '../data/talk_results/pl_T8_no_intervene.txt', '../data/talk_results/powerlawrollout456{0}.txt',
# #                                                 '../data/talk_results/universal_pl_T8_gen4_T4_intervene.txt', '../data/talk_results/powerlawuniversal4_T4{0}_intv.txt',
# #                                                 same_plot=True, normalize_axis_x=False, inset_label=f'$T_g \\rightarrow T_4=0.4$')
# # plt.show()
#
# ## Random Vaccination
# plotting_util.plots_for_nerccs_talk([3, 5, 6, 15], 400,  '../data/talk_results/pl_T8_no_intervene.txt', '../data/talk_results/powerlawrollout456{0}.txt',
#                                                 '../data/talk_results/random_pl_T8_gen4_50pct_intervene.txt', '../data/talk_results/powerlawrandom4{0}_intv.txt',
#                                                 same_plot=True, normalize_axis_x=False, inset_label=f'$T_g \\rightarrow T_4=0$ for \n 50\% Vaccinated')
# plt.show()
#
# ## Random Rollout vaccination
# plotting_util.plots_for_nerccs_talk([3, 10, 15, 20], 400,  '../data/talk_results/pl_T8_no_intervene.txt', '../data/talk_results/powerlawrollout456{0}.txt',
#                                                 '../data/talk_results/rollout_pl_T8_gen456_50pct_intervene.txt', '../data/talk_results/powerlawrollout456{0}_intv.txt',
#                                                 same_plot=True, normalize_axis_x=False, inset_label=f'$T_g \\rightarrow T_4=0$, for 10\% Vaccinated\n $\\rightarrow T_5=0$, for 30\% Vaccinated\n$\\rightarrow T_6=0$, for 50\% Vaccinated')
# plt.show()

## Theoretical Data for the talk
# pgf_formalism.compute_phase_space(21, 400, power_law_dd, 0.8, True, [2, 3, 4, 5, 6, 10, 12, 15, 20], '../data/talk_results/powerlawuniversal4_T4{0}',
#                          intervention_gen=4, intervention_trans=0.4, do_non_interv=False, do_interv=True, intervention_type="universal_intervention")
# plotting_util.outbreak_size_curves([2, 4, 6, 10, 15], 200, '../data/talk_results/powerlawrollout456{0}.txt', '../data/talk_results/powerlawuniversal4_T4{0}_intv.txt', same_plot=True)





