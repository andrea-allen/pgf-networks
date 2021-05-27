#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 10:05:06 2021

@author: nicholasjroberts
"""
import time

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy import stats
from scipy.optimize import root

import src.degree_distributions as degree_distributions
import src.pgf_formalism as pgf_formalism
import src.plotting_util as plotting_util

size = 250

def exp_val_dd(dist):
    ev = 0
    for k in range(len(dist)):
        ev += k*dist[k]
    return ev

power_law_dd = degree_distributions.power_law_degree_distrb(size, mu=5)
extnct_cake_PL_1, psi_1 = pgf_formalism.compute_extinct_prob_all(power_law_dd, T=0.8, n_gens=11, renorm=True)

power_law_dd = degree_distributions.power_law_degree_distrb(size, mu=7.5)
extnct_cake_PL_2, psi_2 = pgf_formalism.compute_extinct_prob_all(power_law_dd, T=0.8, n_gens=11, renorm=True)

power_law_dd = degree_distributions.power_law_degree_distrb(size, mu=10)
extnct_cake_PL_3, psi_3 = pgf_formalism.compute_extinct_prob_all(power_law_dd, T=0.8, n_gens=11, renorm=True)


# generate time series data for pwr law extinction
# =============================================================================
# extnct_cake_PL_1, psi_1 = pgf_formalism.compute_extinct_prob_all(power_law_dd, T=0.8, n_gens=20, renorm=False)
# extnct_ts_1 = np.sum(extnct_cake_PL_1, axis=(1,2))
# =============================================================================

"""
Functions to calculate additional distributions with same mean as power law dd
and their corresponding extinction time series
"""
# =============================================================================
# binom_dd = np.zeros(size)
# p = exp_val_dd(power_law_dd) / size
# for k in range(size):
#     binom_dd[k] = stats.binom.pmf(k, size, p)
# binom_dd = binom_dd / np.sum(binom_dd)
# 
# def find_deg_prob(x, val): return 2*x + (1-x) * 5 - val
# p_small = root(find_deg_prob, 0.5, args=(exp_val_dd(power_law_dd)))['x'][0]
# bi_mod_dd = np.zeros(size)
# bi_mod_dd[2] = p_small
# bi_mod_dd[5] = 1 - p_small
# 
# extnct_cake_Bin_1, psi_2 = pgf_formalism.compute_extinct_prob_all(binom_dd, T=0.8, n_gens =20, renorm=False)
# extnct_ts_2 = np.sum(extnct_cake_Bin_1, axis=(1,2))
# 
# extnct_cake_BiMod_1, psi_3 = pgf_formalism.compute_extinct_prob_all(bi_mod_dd, T=0.8, n_gens =20, renorm=False)
# extnct_ts_3 = np.sum(extnct_cake_BiMod_1, axis=(1,2))
# 
# """
# Plot extinction time series
# """
# fig, ax = plt.subplots()
# ax.scatter(range(20), extnct_ts_1, label = "pwr law")
# ax.scatter(range(20), extnct_ts_2, label = "binom")
# ax.scatter(range(20), extnct_ts_3, label = "bi modal")
# plt.xlabel("generation", fontsize=18)
# plt.ylabel("prob extnct after gen g", fontsize=18)
# plt.legend()
# plt.xticks(range(0,20,2),fontsize=14)
# plt.yticks(fontsize=14)
# plt.tight_layout()
# plt.show()
# =============================================================================


# extnct_cake_PL_2, psi_new = pgf_formalism.compute_extinct_prob_all(power_law_dd, T=0.8, n_gens=11, renorm=True)
#plotting_util.plot_extinct_prob(extnct_cake_Bin_1, x_ax="cumulative infections, s", y_ax="Prob extinction after gen G")
#plotting_util.plot_extinct_prob(extnct_cake_BiMod_1, x_ax="cumulative infections, s", y_ax="Prob extinction after gen G")
plotting_util.plot_extinct_prob(extnct_cake_PL_1, x_ax="cumulative infections, s", y_ax="Prob extinction after gen G", hide_zero=True)
plotting_util.plot_extinct_prob(extnct_cake_PL_2, x_ax="cumulative infections, s", y_ax="Prob extinction after gen G", hide_zero=True)
plotting_util.plot_extinct_prob(extnct_cake_PL_3, x_ax="cumulative infections, s", y_ax="Prob extinction after gen G", hide_zero=True)