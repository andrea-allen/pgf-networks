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

"""
EXAMPLES
###########
#targeted intervention simulations on powerlaw
## SET UP YOUR OWN FILE PATHS BELOW
"""

# CODE TO MODEL
pgf_formalism.compute_phase_space(20, 400, geometric,  0.8, True,
                                  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 17, 18, 19],
                                  f'../nerccs2022/geo_03-24-22_targeted_132pct_{{0}}',
                                  rollout_dict={4:.01, 6:.01, 8:.01, 10:0.01, 12: 0.01}, do_non_interv=False, do_interv=True,
                                  intervention_type="targeted_rollout", pre_vax_correction=True)


