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

import degree_distributions
import pgf_formalism
import plotting_util

size = 400

def exp_val_dd(dist):
    ev = 0
    for k in range(len(dist)):
        ev += k*dist[k]
    return ev

power_law_dd = degree_distributions.power_law_degree_distrb(size)

# =============================================================================
# extnct_cake_PL_1, psi_1 = pgf_formalism.compute_extinct_prob_all(power_law_dd, T=0.8, n_gens=20, renorm=False)
# extnct_ts_1 = np.sum(extnct_cake_PL_1, axis=(1,2))
# 
# binom_dd = np.zeros(size)
# p = exp_val_dd(power_law_dd) / size
# for k in range(size):
#     binom_dd[k] = stats.binom.pmf(k, size, p)
# extnct_cake_B_1, psi_2 = pgf_formalism.compute_extinct_prob_all(binom_dd, T=0.8, n_gens =20, renorm=False)
# extnct_ts_2 = np.sum(extnct_cake_B_1, axis=(1,2))
# 
# fig, ax = plt.subplots()
# ax.scatter(range(20), extnct_ts_1, label = "pwr law")
# ax.scatter(range(20), extnct_ts_2, label = "binom")
# plt.xlabel("generation", fontsize=18)
# plt.ylabel("prob extnct after gen g", fontsize=18)
# plt.legend()
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# plt.tight_layout()
# plt.show()
# 
# =============================================================================

extnct_cake_PL_2, psi_new = pgf_formalism.compute_extinct_prob_all(power_law_dd, T=0.8, n_gens=11, renorm=True)
plotting_util.plot_extinct_prob(extnct_cake_PL_2, x_ax="cumulative infections, s", y_ax="Prob extinction after gen G")