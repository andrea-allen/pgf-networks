"""
File to archive the code used in producing the figures from paper-copies 1
Depends on folder paper1data/ in this same directory.
"""

import matplotlib.pyplot as plt
from src import degree_distributions
from src import plotting_util
from src import figs_for_paper
from analysis import covid
import zipfile

### ARCHIVED PAPER 1 DATA
## Unzip paper1data.zip to run the following figures.
## WARNING: Contents is 2.5 gigabytes

### FIG 1, network cartoon
figs_for_paper.network_drawing()

### FIG 2, cumulative infection distibutions for power law and erdos-renyi networks, respectively:

power_law_q3 = degree_distributions.power_law_degree_distrb(2000, mu=10)  # q is 3, k is 1.7
plotting_util.plot_sims_vs_analytical_multigens([3, 4, 6, 10], 400,
                                                f'paper1data/plaw_T8_10k_120ksims_q3_generational.txt',
                                                'paper1data/powerlaw_q3_T8_ifft_g{0}.txt',
                                                inset_to_plot=power_law_q3, inset_title='$p_k = k^{-2}e^{-k/10}$',
                                                same_plot=True, normalize_axis_x=False, plot_distribution_inset=True)
plt.show()
# # #
poisson = degree_distributions.binomial_degree_distb(2000, 2.5)
plotting_util.plot_sims_vs_analytical_multigens([3, 4, 6, 10], 400,
                                                f'./paper1data/erdos_renyi_10k_60ksims_combo_generational.txt',
                                                './paper1data/ER_q2.5_T8_ifft_g{0}.txt', inset_to_plot=poisson,
                                                inset_title='$p_k \\approx \\frac{\\lambda^ke^{-\\lambda}}{k!}, \\lambda=2.5$',
                                                same_plot=True, normalize_axis_x=False, plot_distribution_inset=True)
plt.show()

### FIG 3
figs_for_paper.results_plots('./paper1data/plaw_06_10_40knet', q_degree=3.04, active_gen_sizes_on=True)

### FIG 4
covid.contour_figs('./paper1data')
plt.show()

### FIG 5
covid.contour_fig3()
plt.show()