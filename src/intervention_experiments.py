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

LOCAL_DATA_PATH = "../data"

power_law_q2 = degree_distributions.power_law_degree_distrb(400, mu=10)
k_mean_degree = 2.5
# k_mean_degree = 5
er_degree_dist = degree_distributions.binomial_degree_distb(400, k_mean_degree)

def plot_manatee_concept():
    ## Here: 2/8/22: Plotting just the transmission expression to understand if it's what we want:
    plt.figure('one intervention discrepancy')
    beta = 0.8
    gamma = 0.2
    # T0 = beta/(beta+gamma)
    betas = np.full(40, beta)
    gammas = np.full(40, gamma)
    qs = np.full(40, 2.5)
    vacc = [0 for i in range(40)]
    vacc[3:] = [0.5 for i in range(40-3)]
    t_values = np.zeros(40)
    for g in range(40):
        t_values[g] = mnte.t_of_g(betas, gammas, qs, vacc, g)

    plt.plot(np.arange(40), t_values, label='40 gens in the series')
    plt.text(3, (.5*beta/(beta+gamma)), 'intervention')
    plt.scatter([3], [(.5*beta/(beta+gamma))])

    betas = np.full(20, beta)
    gammas = np.full(20, gamma)
    qs = np.full(20, 2.5)
    vacc = [0 for i in range(20)]
    vacc[3:] = [0.5 for i in range(20-3)]
    t_values = np.zeros(20)
    for g in range(20):
        t_values[g] = mnte.t_of_g(betas, gammas, qs, vacc, g)

    plt.plot(np.arange(20), t_values, label='20 gens in the series')
    plt.text(3, (.5*beta/(beta+gamma)), 'intervention')
    plt.scatter([3], [(.5*beta/(beta+gamma))])
    plt.legend()
    plt.ylabel('Transmissibility T(g)')
    plt.xlabel('Generation g')
    plt.plot(np.arange(3, 40), np.full(37,(.5*beta/(beta+gamma))), label='T(g) should be', ls='--', color='k')
    # plt.show()

    # version that's no intervention
    plt.figure('no intervention')
    betas = np.full(50, beta)
    gammas = np.full(50, gamma)
    qs = np.full(50, 2.5)
    vacc = [0 for i in range(50)]
    # vacc[3:] = [0.5 for i in range(50-3)]
    t_values = np.zeros(50)
    for g in range(50):
        t_values[g] = mnte.t_of_g(betas, gammas, qs, vacc, g)

    plt.plot(np.arange(50), t_values, label='baseline, no intervention, 50 gens')

    # AGAIN, with a longer power series:
    betas = np.full(20, beta)
    gammas = np.full(20, gamma)
    qs = np.full(20, 2.5)
    vacc = [0 for i in range(20)]
    # vacc[3:] = [0.5 for i in range(20-3)]
    t_values = np.zeros(20)
    for g in range(20):
        t_values[g] = mnte.t_of_g(betas, gammas, qs, vacc, g)

    plt.plot(np.arange(20), t_values, label='baseline, no intervention, 20 gens')
    plt.plot(np.arange(0, 50), np.full(50,(beta/(beta+gamma))), label='T(g) should be', ls='--', color='k')
    # plt.text(3, (beta/(beta+gamma)), 'intervention')
    # plt.scatter([3], [(beta/(beta+gamma))])
    plt.ylabel('Transmissibility T(g)')
    plt.xlabel('Generation g')
    plt.legend()

    plt.figure('multiple interventions')
    betas = np.full(200, beta)
    gammas = np.full(200, gamma)
    qs = np.full(200, 2.5)
    vacc = [0 for i in range(200)]
    vacc[3:] = [0.2 for i in range(200-3)]
    vacc[7:] = [0.5 for i in range(200-7)]
    t_values = np.zeros(200)
    for g in range(200):
        t_values[g] = mnte.t_of_g(betas, gammas, qs, vacc, g)

    plt.plot(np.arange(200), t_values, label='2 interventions, 200 gens')

    # AGAIN but with 20 gens
    betas = np.full(20, beta)
    gammas = np.full(20, gamma)
    qs = np.full(20, 2.5)
    vacc = [0 for i in range(20)]
    vacc[3:] = [0.2 for i in range(20-3)]
    vacc[7:] = [0.5 for i in range(20-7)]
    t_values = np.zeros(20)
    for g in range(20):
        t_values[g] = mnte.t_of_g(betas, gammas, qs, vacc, g)

    plt.plot(np.arange(20), t_values, label='2 interventions, 20 gens')
    plt.plot(np.arange(0, 200), np.full(200,(.5*beta/(beta+gamma))), label='T(g) should be', ls='--', color='k')

    # plt.text(3, (beta/(beta+gamma)), 'intervention')
    # plt.scatter([3], [(beta/(beta+gamma))])
    plt.ylabel('Transmissibility T(g)')
    plt.xlabel('Generation g')
    plt.legend()

    plt.show()

############
# new example with a slower rollout for random vaccination, with powerlaw
# ensemble.run_ensemble_intervention_effects(power_law_q2, '../data/rollouts/random/pl_02-22-22_g5-9-13', 6000, 10000,
#                                            init_T=0.8, intervention_gen_list=[5, 9, 13], prop_reduced_list=[.2, .2, .1],
#                                            intervention_type="random_rollout", run_regular=False)

# pgf_formalism.compute_phase_space(20, 400, power_law_q2,  0.8, True,
#                                   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 17, 18, 19],
#                                   f'{LOCAL_DATA_PATH}/rollouts/random/pl_02-22-22_g5-9-13_{{0}}',
#                                   rollout_dict={5: .2, 9: .2, 13: .1}, do_non_interv=False, do_interv=True,
#                                   intervention_type="random_rollout", pre_vax_correction=True)

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
plt.show()

############



# ensemble.run_ensemble_intervention_effects(er_degree_dist, '../data/rollouts/random/er_g5', 20000, 10000,
#                                            init_T=0.8, intervention_gen_list=[3, 4, 6], prop_reduced_list=[.2, .2, .3],
#                                            intervention_type="random_rollout", run_regular=False)

# ER with the magic random rollout effects
pgf_formalism.compute_phase_space(20, 400, er_degree_dist, 0.8, True,
                                  [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 16, 17, 18, 19],
                                  f'{LOCAL_DATA_PATH}/rollouts/random/er_manatee_02_15_alt3_{{0}}',
                                  rollout_dict={3: .2, 4: .2, 6: .3}, do_non_interv=False, do_interv=True,
                                  intervention_type="random_rollout", pre_vax_correction=True)

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