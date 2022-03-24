
from src import degree_distributions

from src import pgf_formalism


LOCAL_DATA_PATH = "../data"

power_law_q2 = degree_distributions.power_law_degree_distrb(400, mu=10)
k_mean_degree = 2.5
# k_mean_degree = 5
er_degree_dist = degree_distributions.binomial_degree_distb(400, k_mean_degree)

geom_degree_dist = degree_distributions.exponential_degree_dist(400, .6) # n and then p
"""
EXAMPLES
###########
#targeted intervention simulations on powerlaw
## SET UP YOUR OWN FILE PATHS BELOW
"""

# CODE TO MODEL
pgf_formalism.compute_phase_space(20, 400, geom_degree_dist,  0.8, True,
                                  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19],
                                  f'../nerccs2022/geo_03-24-22_targeted_slow_half_{{0}}',
                                  rollout_dict={4:.1, 6:.1, 8:.1, 10:.1, 12: .1}, do_non_interv=False, do_interv=True,
                                  intervention_type="targeted_rollout", pre_vax_correction=True)


