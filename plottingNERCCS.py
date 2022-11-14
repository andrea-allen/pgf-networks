#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 15:04:49 2022

@author: mariahboudreau
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300

plt.rcParams["font.family"] = "Times New Roman"
x_start = 1
x_lim = 400

### Plotting for NERCCS

fig, ax1 = plt.subplots()
ax1.semilogy()
ax1.set_ylim(.00005, .1)  
ax1.set_xlabel("Cumulative infections $s$ ")
ax1.set_ylabel("Probability")

# data = np.loadtxt("random/sims_03-21-22_rand20k_intervene.txt", delimiter=',')
# count = 0
# for gen in [2,4,6,8]:
#     count += 1
#     time_series = data[gen][x_start:x_lim]
#     for x in range(0, x_lim-x_start-1):
#         if time_series[x]<10**(-5):
#             print(time_series[x])
#             time_series[x] = (time_series[x-1] + time_series[x+1])/2
      
#     ax1.plot(x_vals, time_series, label= f"gen {gen}")


import seaborn as sns
color = sns.color_palette("mako", n_colors=10)
        
#### Phase space       

fig, (ax1, ax2, ax3) = plt.subplots(3,1, sharey= True, sharex = True)
fig.set_size_inches(5, 15)
ax1.semilogy()
ax1.set_ylim(.00005, .1)  
ax1.set_ylabel("Probability")
ax2.set_ylabel("Probability")
ax3.set_ylabel("Probability")
ax3.set_xlabel("Cumulative infections, $s$")


 



for i, c in zip([2,4,6,8], color):
    psi_g = np.loadtxt(f"nerccs2022/geo_05-5-22_targeted_rollout_2.5_resultsfigs_{i}_intv.txt", delimiter=',') #Make sure that this is the right relative path
    inverted_s_m = psi_g.T
    ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
    ps_g_analytical = ps_g_analytical / np.sum(ps_g_analytical)  # normalize
    x_vals = np.arange(1, 400)
    ax3.plot(x_vals, ps_g_analytical[1:400], c = c)
    

for i, c in zip([2,4,6,8], color):
    psi_g = np.loadtxt(f"nerccs2022/geo_05-5-22_random_rollout_15_{i}_intv.txt", delimiter=',')
    inverted_s_m = psi_g.T
    ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
    ps_g_analytical = ps_g_analytical / np.sum(ps_g_analytical)  # normalize
    x_vals = np.arange(1, 400)
    ax3.plot(x_vals, ps_g_analytical[1:400], c = c, linestyle = "dashdot")
      
    
# ax3.semilogy()
# ax3.set_ylim(.00005, .1) 

# for i, c in zip([2,4,6,8,10], color):
#     psi_g = np.loadtxt(f"nerccs2022/geo_05-5-22_targeted_rollout_5_resultsfigs_{i}_intv.txt", delimiter=',') #Make sure that this is the right relative path
#     inverted_s_m = psi_g.T
#     ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
#     ps_g_analytical = ps_g_analytical / np.sum(ps_g_analytical)  # normalize
#     x_vals = np.arange(1, 400)
#     ax3.plot(x_vals, ps_g_analytical[1:400], c = c)
      
# for i, c in zip([2,4,6,8,10], color):
#     psi_g = np.loadtxt(f"nerccs2022/geo_05-5-22_random_rollout_15_{i}_intv.txt", delimiter=',')
#     inverted_s_m = psi_g.T
#     ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
#     ps_g_analytical = ps_g_analytical / np.sum(ps_g_analytical)  # normalize
#     x_vals = np.arange(1, 400)
#     ax3.plot(x_vals, ps_g_analytical[1:400], c = c, linestyle = "dashdot")
    
ax2.semilogy()
ax2.set_ylim(.00005, .1) 

for i, c in zip([2,4,6,8], color):
    psi_g = np.loadtxt(f"nerccs2022/geo_05-5-22_targeted_rollout_7.5_resultsfigs_{i}_intv.txt", delimiter=',') #Make sure that this is the right relative path
    inverted_s_m = psi_g.T
    ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
    ps_g_analytical = ps_g_analytical / np.sum(ps_g_analytical)  # normalize
    x_vals = np.arange(1, 400)
    ax2.plot(x_vals, ps_g_analytical[1:400], c = c)
      
for i, c in zip([2,4,6,8,10], color):
    psi_g = np.loadtxt(f"nerccs2022/geo_05-5-22_random_rollout_15_{i}_intv.txt", delimiter=',')
    inverted_s_m = psi_g.T
    ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
    ps_g_analytical = ps_g_analytical / np.sum(ps_g_analytical)  # normalize
    x_vals = np.arange(1, 400)
    ax2.plot(x_vals, ps_g_analytical[1:400], c = c, linestyle = "dashdot")
    
    
ax1.semilogy()
ax1.set_ylim(.00005, .1) 

for i, c in zip([2,4,6,8], color):
    psi_g = np.loadtxt(f"nerccs2022/geo_05-5-22_targeted_rollout_10_resultsfigs_{i}_intv.txt", delimiter=',') #Make sure that this is the right relative path
    inverted_s_m = psi_g.T
    ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
    ps_g_analytical = ps_g_analytical / np.sum(ps_g_analytical)  # normalize
    x_vals = np.arange(1, 400)
    ax1.plot(x_vals, ps_g_analytical[1:400], c = c, label = f"Gen {i} (targeted)")
      


for i, c in zip([2,4,6,8], color):
    psi_g = np.loadtxt(f"nerccs2022/geo_05-5-22_random_rollout_15_{i}_intv.txt", delimiter=',')
    inverted_s_m = psi_g.T
    ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
    ps_g_analytical = ps_g_analytical / np.sum(ps_g_analytical)  # normalize
    x_vals = np.arange(1, 400)
    ax1.plot(x_vals, ps_g_analytical[1:400], c = c, linestyle = "dashdot", label = f"Gen {i} (random)")

fig.legend(loc = 'upper center' , frameon=False, ncol = 2, fontsize = 12)




#%%%% Plotting individual (random with sims at gen 5)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns

plt.rcParams["font.family"] = "Times New Roman"
x_start = 1
x_lim = 400
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
fig, ax1 = plt.subplots()
ax1.semilogy()
ax1.set_ylim(.00005, .1)  
ax1.set_xlabel("Cumulative infections, $s$ ")
ax1.set_ylabel("Probability")
color = sns.color_palette("mako", n_colors=6)

data = np.loadtxt("/geo_09-21-22_random_gen5_sims_intervene.txt", delimiter=',')
x_vals = np.arange(1, 400)

for gen, c in zip([2,4,6,8,10,12,14,16], color):
    
    time_series = data[gen][x_start:x_lim]
    for x in range(0, x_lim-x_start-1):
        if time_series[x]<10**(-5):
            print(time_series[x])
            time_series[x] = (time_series[x-1] + time_series[x+1])/2
    
    if gen == 10:
        print(time_series)
    ax1.plot(x_vals, time_series, c = c, label= f"Gen {gen}")




for i, c in zip([2,4,6,8,10,12,14,16], color):
    psi_g = np.loadtxt(f"geo_09-21-22_random_gen5_25_{i}_intv.txt", delimiter=',')
    inverted_s_m = psi_g.T
    ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
    ps_g_analytical = ps_g_analytical / np.sum(ps_g_analytical)  # normalize

    x_vals = np.arange(1, 400)
    ax1.plot(x_vals, ps_g_analytical[1:400], c = c)

plt.legend(frameon = False, fontsize = 8)

plt.show()


#%%%% Plotting individual (targeted with sims at gen 5) 


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
plt.rcParams["font.family"] = "Times New Roman"
x_start = 1
x_lim = 400

import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300

fig, ax1 = plt.subplots()
ax1.semilogy()
ax1.set_ylim(.00005, .1)  
ax1.set_xlabel("Cumulative infections $s$ ")
ax1.set_ylabel("Probability")
color = sns.color_palette("mako", n_colors=6)

data = np.loadtxt("nerccs2022/geo_10-10-22_targeted_gen5_sims_p015_intervene.txt", delimiter=',')
x_vals = np.arange(1, 400)

for gen, c in zip([2,4,6,8,10,12,14], color):
    
    time_series = data[gen][x_start:x_lim]
    for x in range(0, x_lim-x_start-1):
        if time_series[x]<10**(-5):
            print(time_series[x])
            time_series[x] = (time_series[x-1] + time_series[x+1])/2
      
    ax1.plot(x_vals, time_series, c = c, label= f"Gen {gen}")




for i, c in zip([2,4,6,8,10,12,14], color):
    psi_g = np.loadtxt(f"nerccs2022/geo_10-5-22_targeted_gen5_p015_{i}_intv.txt", delimiter=',')
    inverted_s_m = psi_g.T
    ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
    ps_g_analytical = ps_g_analytical / np.sum(ps_g_analytical)  # normalize
    x_vals = np.arange(1, 400)
    ax1.plot(x_vals, ps_g_analytical[1:400], c = c)

plt.legend(frameon = False, fontsize = 8)



#%%%% Plotting individual (random rollout with sims at 0.25%)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
plt.rcParams["font.family"] = "Times New Roman"
x_start = 1
x_lim = 400

import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300

fig, ax1 = plt.subplots()
ax1.semilogy()
ax1.set_ylim(.00005, .1)  
ax1.set_xlabel("Cumulative infections $s$ ")
ax1.set_ylabel("Probability")
color = sns.color_palette("mako", n_colors=6)

data = np.loadtxt("nerccs2022/geo_05-5-22_random_rollout_005_sims_intervene.txt", delimiter=',')
x_vals = np.arange(1, 400)

for gen, c in zip([2,4,6,8,10,12,14,16], color):
    
    time_series = data[gen][x_start:x_lim]
    for x in range(0, x_lim-x_start-1):
        if time_series[x]<10**(-5):
            print(time_series[x])
            time_series[x] = (time_series[x-1] + time_series[x+1])/2
      
    ax1.plot(x_vals, time_series, c = c, label= f"Gen {gen}")




for i, c in zip([2,4,6,8,10,12,14,16], color):
    psi_g = np.loadtxt(f"nerccs2022/geo_05-11-22_random_rollout_005_{i}_intv.txt", delimiter=',')
    inverted_s_m = psi_g.T
    ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
    ps_g_analytical = ps_g_analytical / np.sum(ps_g_analytical)  # normalize
    x_vals = np.arange(1, 400)
    ax1.plot(x_vals, ps_g_analytical[1:400], c = c)

plt.legend(frameon = False, fontsize = 8)





#%%%% Plotting individual (targeted with sims rollout cumulative 2%)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
plt.rcParams["font.family"] = "Times New Roman"
x_start = 1
x_lim = 400

import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300

fig, ax1 = plt.subplots()
ax1.semilogy()
ax1.set_ylim(.00005, .1)  
ax1.set_xlabel("Cumulative infections $s$ ")
ax1.set_ylabel("Probability")
color = sns.color_palette("mako", n_colors=6)

data = np.loadtxt("nerccs2022/geo_09-27-22_75k_targeted_rollout_2_sims_intervene.txt", delimiter=',') # initial targeted at .15%
x_vals = np.arange(1, 400)

for gen, c in zip([2,4,6,8,10,12], color):
    
    time_series = data[gen][x_start:x_lim]
    for x in range(0, x_lim-x_start-1):
        if time_series[x]<10**(-5):
            print(time_series[x])
            time_series[x] = (time_series[x-1] + time_series[x+1])/2
      
    ax1.plot(x_vals, time_series, c = c, label= f"Gen {gen}")




for i, c in zip([2,4,6,8,10,12], color):
    psi_g = np.loadtxt(f"nerccs2022/geo_10-1-22_targeted_rollout_2_{i}_intv.txt", delimiter=',')
    inverted_s_m = psi_g.T
    ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
    ps_g_analytical = ps_g_analytical / np.sum(ps_g_analytical)  # normalize
    x_vals = np.arange(1, 400)
    ax1.plot(x_vals, ps_g_analytical[1:400], c = c)

plt.legend(frameon = False, fontsize = 8)


#%%% Setting up the functions for metrics

import numpy as np
import matplotlib.pyplot as plt

#In case we get issues with finding exact values 
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

#In case we get issues with finding exact values 
def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

### METRICS

# Average cases 
def average_cases(dist):
    temp = 0
    for i in range(len(dist)):
        temp += dist[i]*i
        
    return temp


#  Worst Case scenerio
# Max cases such that P(c) is greater than a threshold
def worst_case_scenario(dist, thresh):
    worst_case = find_nearest_idx(dist,thresh)

    return worst_case

# Shit hits the fan (what is the prob of this many cases happening)
def s_hits_fan(dist, s_cases):
    
    return dist[s_cases]

# Probability of more cases with intervention
def pointless_interv(interv_dist, non_dist, max_cases):
    #Probability of cases from targeted are higher than non intervention
    # this is calculated by the sum of Prob of cases( targeted) * CDF(cases(none)) for all cases
    pointless_metric = np.zeros((max_cases))
    
    #CDF calculation
    cdf = 1. * np.arange(len(non_dist))/(len(non_dist))
    
    for i in range(max_cases):
        pointless_metric[i] = interv_dist[i]*cdf[i]
        
        
    return np.sum(pointless_metric)



#%%% Plotting the comparison metrics
    
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300

plt.rcParams["font.family"] = "Times New Roman"
fig, axs = plt.subplots(4,2, sharex='col', sharey='row')

fig.set_size_inches(9, 12)

case_vals = [0, .025, .05, .075, .10, .125, .15]

# 0.001, 0.0025, 0.005, 0.0075, 0.010, 0.0125, 0.015, 0.0175, 0.020

points = len(case_vals)

gen5_target_per = np.zeros((points))

gen10_target_per = np.zeros((points))

avg_cases_gen5 = np.zeros((points))

avg_cases_gen10 = np.zeros((points))

worst_gen5 = np.zeros((points))

worst_gen10 = np.zeros((points))

fan_gen5 = np.zeros((points))

fan_gen10 = np.zeros((points))

pointless_gen5 = np.zeros((points))

pointless_gen10 = np.zeros((points))

thresh_worst = 10**(-3)

thresh_fan = 500

thresh_point = 799

# Keep this consistant for generation 5 and generation 10

psi_g = np.loadtxt(f"nerccs2022/geo_05-11-22_random_rollout_005_5_intv.txt", delimiter=',')
inverted_s_m = psi_g.T
ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
random_constant_gen5 = ps_g_analytical / np.sum(ps_g_analytical)  # normalize

psi_g = np.loadtxt(f"nerccs2022/geo_05-11-22_random_rollout_005_10_intv.txt", delimiter=',')
inverted_s_m = psi_g.T
ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
random_constant_gen10 = ps_g_analytical / np.sum(ps_g_analytical)  # normalize

# Non intervention

psi_g = np.loadtxt(f"nerccs2022/geo_10-3-22_non_interv_5.txt", delimiter=',') 
inverted_s_m = psi_g.T
ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
gen5_none = ps_g_analytical / np.sum(ps_g_analytical)  # normalize

avg_cases_gen5[0] = average_cases(gen5_none)
worst_gen5[0] = worst_case_scenario(gen5_none, thresh_worst)
fan_gen5[0] = s_hits_fan(gen5_none, thresh_fan)
pointless_gen5[0] = pointless_interv(gen5_none, gen5_none, thresh_point)

psi_g = np.loadtxt(f"nerccs2022/geo_10-3-22_non_interv_10.txt", delimiter=',') 
inverted_s_m = psi_g.T
ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
gen10_none = ps_g_analytical / np.sum(ps_g_analytical)  # normalize

avg_cases_gen10[0] = average_cases(gen10_none)
worst_gen10[0] = worst_case_scenario(gen10_none, thresh_worst)
fan_gen10[0] = s_hits_fan(gen10_none, thresh_fan)
pointless_gen10[0] = pointless_interv(gen10_none, gen10_none, thresh_point)



for i, percent in zip(range(1,points),case_vals[1:]):
    
    psi_g = np.loadtxt(f"nerccs2022/geo_10-3-22_targeted_rollout_{percent}_resultsfigs_5_intv.txt", delimiter=',') #Make sure that this is the right relative path
    inverted_s_m = psi_g.T
    ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
    gen5_target_per = np.array(ps_g_analytical / np.sum(ps_g_analytical))  
    
    
    avg_cases_gen5[i] = average_cases(gen5_target_per)
    worst_gen5[i] = worst_case_scenario(gen5_target_per, thresh_worst)
    fan_gen5[i] = s_hits_fan(gen5_target_per, thresh_fan)
    pointless_gen5[i] = pointless_interv(gen5_target_per, gen5_none, thresh_point)



for i, percent in zip(range(1, points),case_vals[1:]):
    
    psi_g = np.loadtxt(f"nerccs2022/geo_10-3-22_targeted_rollout_{percent}_resultsfigs_10_intv.txt", delimiter=',') #Make sure that this is the right relative path
    inverted_s_m = psi_g.T
    ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
    
    gen10_target_per = np.array(ps_g_analytical / np.sum(ps_g_analytical))  
    avg_cases_gen10[i] = average_cases(gen10_target_per)
    worst_gen10[i] = worst_case_scenario(gen10_target_per, thresh_worst)
    fan_gen10[i] = s_hits_fan(gen10_target_per, thresh_fan)
    pointless_gen10[i] = pointless_interv(gen10_target_per, gen10_none, thresh_point)
    




rand_gen5_line = average_cases(random_constant_gen5)
rand_gen10_line = average_cases(random_constant_gen10)

avg_bars_gen5 = rand_gen5_line - avg_cases_gen5
bar_colors_gen5 = []
bar_start_gen5 = np.zeros((len(avg_bars_gen5)))
bar_end_gen5 = np.zeros((len(avg_bars_gen5)))



for num, c in zip(avg_bars_gen5, range(len(avg_bars_gen5))):
    if num < 0:
        bar_colors_gen5.append("red")
        bar_start_gen5[c] = avg_cases_gen5[c]
        bar_end_gen5[c] = rand_gen5_line
    else:
        bar_colors_gen5.append("green")
        bar_start_gen5[c] = rand_gen5_line
        bar_end_gen5[c] =  avg_cases_gen5[c]
        
avg_bars_gen10 = rand_gen10_line - avg_cases_gen10
bar_colors_gen10 = []
bar_start_gen10 = np.zeros((len(avg_bars_gen5)))
bar_end_gen10 = np.zeros((len(avg_bars_gen5)))

for num, c in zip(avg_bars_gen10, range(len(avg_bars_gen10))):
    if num < 0:
        bar_colors_gen10.append("red")
        bar_start_gen10[c] = avg_cases_gen10[c]
        bar_end_gen10[c] = rand_gen10_line
    else:
        bar_colors_gen10.append("green")
        bar_start_gen10[c] = rand_gen10_line
        bar_end_gen10[c] =  avg_cases_gen10[c]

axs[0,0].set_ylabel("Number of cases")
axs[0,0].set_title("Average cases at Gen. 5")
axs[0,1].set_title("Average cases at Gen. 10")

axs[0,0].vlines(case_vals, bar_end_gen5, bar_start_gen5, colors = bar_colors_gen5)
axs[0,1].vlines(case_vals, bar_end_gen10, bar_start_gen10, colors = bar_colors_gen10)

axs[0,0].scatter(case_vals, avg_cases_gen5, color = bar_colors_gen5)
axs[0,0].axhline(rand_gen5_line, linestyle = 'dashed', color = "black")

axs[0,1].scatter(case_vals, avg_cases_gen10, color = bar_colors_gen5)
axs[0,1].axhline(rand_gen10_line, linestyle = 'dashed', color = "black")



rand_gen5_line = worst_case_scenario(random_constant_gen5, thresh_worst)
rand_gen10_line = worst_case_scenario(random_constant_gen10, thresh_worst)


worst_bars_gen5 = rand_gen5_line - worst_gen5
bar_colors_gen5 = []
bar_start_gen5 = np.zeros((len(worst_bars_gen5)))
bar_end_gen5 = np.zeros((len(worst_bars_gen5)))

for num, c in zip(worst_bars_gen5, range(len(worst_bars_gen5))):
    if num < 0:
        bar_colors_gen5.append("red")
        bar_start_gen5[c] = worst_gen5[c]
        bar_end_gen5[c] = rand_gen5_line
    else:
        bar_colors_gen5.append("green")
        bar_start_gen5[c] = rand_gen5_line
        bar_end_gen5[c] =  worst_gen5[c]
        
worst_bars_gen10 = rand_gen10_line - worst_gen10
bar_colors_gen10 = []
bar_start_gen10 = np.zeros((len(worst_gen10)))
bar_end_gen10 = np.zeros((len(worst_gen10)))

for num, c in zip(worst_bars_gen10, range(len(worst_bars_gen10))):
    if num < 0:
        bar_colors_gen10.append("red")
        bar_start_gen10[c] = worst_gen10[c]
        bar_end_gen10[c] = rand_gen10_line
    else:
        bar_colors_gen10.append("green")
        bar_start_gen10[c] = rand_gen10_line
        bar_end_gen10[c] =  worst_gen10[c]

axs[1,0].set_ylabel("Number of cases")
axs[1,0].set_title("Best - Worst case at Gen. 5")
axs[1,1].set_title("Best - Worst case at Gen. 10")

axs[1,0].vlines(case_vals, bar_end_gen5, bar_start_gen5, colors = bar_colors_gen5)
axs[1,1].vlines(case_vals, bar_end_gen10, bar_start_gen10, colors = bar_colors_gen10)

axs[1,0].scatter(case_vals, worst_gen5, marker = "^", color = bar_colors_gen5)
axs[1,0].axhline(rand_gen5_line, linestyle = 'dashed', color ="black")

axs[1,1].scatter(case_vals, worst_gen10, marker = "^",color = bar_colors_gen10)
axs[1,1].axhline(rand_gen10_line, linestyle = 'dashed', color = "black")



rand_gen5_line = s_hits_fan(random_constant_gen5, thresh_fan)
rand_gen10_line = s_hits_fan(random_constant_gen10, thresh_fan)


fan_bars_gen5 = rand_gen5_line - fan_gen5
bar_colors_gen5 = []
bar_start_gen5 = np.zeros((len(fan_bars_gen5)))
bar_end_gen5 = np.zeros((len(fan_bars_gen5)))

for num, c in zip(fan_bars_gen5, range(len(fan_bars_gen5))):
    if num < 0:
        bar_colors_gen5.append("red")
        bar_start_gen5[c] = fan_gen5[c]
        bar_end_gen5[c] = rand_gen5_line
    else:
        bar_colors_gen5.append("green")
        bar_start_gen5[c] = rand_gen5_line
        bar_end_gen5[c] =  fan_gen5[c]
        
fan_bars_gen10 = rand_gen10_line - fan_gen10
bar_colors_gen10 = []
bar_start_gen10 = np.zeros((len(fan_bars_gen10)))
bar_end_gen10 = np.zeros((len(fan_bars_gen10)))

for num, c in zip(fan_bars_gen10, range(len(fan_bars_gen10))):
    if num < 0:
        bar_colors_gen10.append("red")
        bar_start_gen10[c] = fan_gen10[c]
        bar_end_gen10[c] = rand_gen10_line
    else:
        bar_colors_gen10.append("green")
        bar_start_gen10[c] = rand_gen10_line
        bar_end_gen10[c] =  fan_gen10[c]

axs[2,0].set_ylabel("Probability")
axs[2,0].semilogy()
axs[2,0].set_title("Critical level of cases at Gen. 5")
axs[2,1].set_title("Critical level of cases at Gen. 10")

axs[2,0].vlines(case_vals, bar_end_gen5, bar_start_gen5, colors = bar_colors_gen5)
axs[2,1].vlines(case_vals, bar_end_gen10, bar_start_gen10, colors = bar_colors_gen10)

axs[2,0].scatter(case_vals, fan_gen5, marker = "s", color = bar_colors_gen5)
axs[2,0].axhline(rand_gen5_line, linestyle = 'dashed', color = "black")

axs[2,1].scatter(case_vals, fan_gen10, marker = "s", color = bar_colors_gen10)
axs[2,1].axhline(rand_gen10_line, linestyle = 'dashed', color = "black")



rand_gen5_line = pointless_interv(random_constant_gen5, gen5_none, thresh_point)
rand_gen10_line = pointless_interv(random_constant_gen10, gen10_none, thresh_point)

point_bars_gen5 = rand_gen5_line - pointless_gen5
bar_colors_gen5 = []
bar_start_gen5 = np.zeros((len(point_bars_gen5)))
bar_end_gen5 = np.zeros((len(point_bars_gen5)))

for num, c in zip(point_bars_gen5, range(len(point_bars_gen5))):
    if num < 0:
        bar_colors_gen5.append("red")
        bar_start_gen5[c] = pointless_gen5[c]
        bar_end_gen5[c] = rand_gen5_line
    else:
        bar_colors_gen5.append("green")
        bar_start_gen5[c] = rand_gen5_line
        bar_end_gen5[c] =  pointless_gen5[c]
        
point_bars_gen10 = rand_gen10_line - pointless_gen10
bar_colors_gen10 = []
bar_start_gen10 = np.zeros((len(point_bars_gen10 )))
bar_end_gen10 = np.zeros((len(point_bars_gen10 )))

for num, c in zip(point_bars_gen10, range(len(point_bars_gen10))):
    if num < 0:
        bar_colors_gen10.append("red")
        bar_start_gen10[c] = pointless_gen10[c]
        bar_end_gen10[c] = rand_gen10_line
    else:
        bar_colors_gen10.append("green")
        bar_start_gen10[c] = rand_gen10_line
        bar_end_gen10[c] =  pointless_gen10[c]



axs[3,0].set_ylabel("Probability")
axs[3,0].set_title("Minimal effect intervention at Gen. 5")
axs[3,1].set_title("Minimal effect intervention at Gen. 10")

axs[3,0].semilogy()

axs[3,0].set_xlabel("% Vaccinated for each rollout")
axs[3,1].set_xlabel("% Vaccinated for each rollout")
axs[3,0].set_xticks(case_vals)
axs[3,1].set_xticks(case_vals)

axs[3,0].vlines(case_vals, bar_end_gen5, bar_start_gen5, colors = bar_colors_gen5)
axs[3,1].vlines(case_vals, bar_end_gen10, bar_start_gen10, colors = bar_colors_gen10)

axs[3,0].scatter(case_vals, pointless_gen5, marker = "*", color = bar_colors_gen5)
axs[3,0].axhline(rand_gen5_line, linestyle = 'dashed', color = "black")

axs[3,1].scatter(case_vals, pointless_gen10, marker = "*", color = bar_colors_gen10)
axs[3,1].axhline(rand_gen10_line, linestyle = 'dashed', color = "black")


#%%% Old plotting

# ######  
# # 2.5% plotting
# psi_g = np.loadtxt(f"pgf-networks/nerccs2022/geo_04-25-22_targeted_rollout_2point5per_5_intv.txt", delimiter=',') #Make sure that this is the right relative path
# inverted_s_m = psi_g.T
# ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
# gen5_target_2point5per = ps_g_analytical / np.sum(ps_g_analytical)  # normalize

# psi_g = np.loadtxt(f"pgf-networks/nerccs2022/geo_04-25-22_targeted_rollout_2point5per_10_intv.txt", delimiter=',') #Make sure that this is the right relative path
# inverted_s_m = psi_g.T
# ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
# gen10_target_2point5per = ps_g_analytical / np.sum(ps_g_analytical)  # normalize

#####
# 5% plotting
psi_g = np.loadtxt(f"pgf-networks/nerccs2022/geo_04-25-22_targeted_rollout_5per_5_intv.txt", delimiter=',') #Make sure that this is the right relative path
inverted_s_m = psi_g.T
ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
gen5_target_5per = ps_g_analytical / np.sum(ps_g_analytical)  # normalize
    
    
psi_g = np.loadtxt(f"pgf-networks/nerccs2022/geo_04-25-22_targeted_rollout_5per_10_intv.txt", delimiter=',') #Make sure that this is the right relative path
inverted_s_m = psi_g.T
ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
gen10_target_5per = ps_g_analytical / np.sum(ps_g_analytical)  # normalize


avg_cases_gen5[1] = average_cases(gen5_target_5per)

avg_cases_gen10[1] = average_cases(gen10_target_5per)






#######
# 10% plotting
psi_g = np.loadtxt(f"pgf-networks/nerccs2022/geo_04-25-22_targeted_rollout_10per_5_intv.txt", delimiter=',') 

inverted_s_m = psi_g.T
ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
gen5_target_10per = ps_g_analytical / np.sum(ps_g_analytical)  # normalize\

psi_g = np.loadtxt(f"pgf-networks/nerccs2022/geo_04-25-22_targeted_rollout_10per_10_intv.txt", delimiter=',') 

inverted_s_m = psi_g.T
ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
gen10_target_10per = ps_g_analytical / np.sum(ps_g_analytical)  # normalize

avg_cases_gen5[2] = average_cases(gen5_target_10per)

avg_cases_gen10[2] = average_cases(gen10_target_10per)



#####
# 20% vaccination
psi_g = np.loadtxt(f"pgf-networks/nerccs2022/geo_05-5-22_targeted_rollout_20per_5_intv.txt", delimiter=',') #Make sure that this is the right relative path
inverted_s_m = psi_g.T
ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
gen5_target_20per = ps_g_analytical / np.sum(ps_g_analytical)  # normalize

psi_g = np.loadtxt(f"pgf-networks/nerccs2022/geo_05-5-22_targeted_rollout_20per_10_intv.txt", delimiter=',') #Make sure that this is the right relative path
inverted_s_m = psi_g.T
ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
gen10_target_20per = ps_g_analytical / np.sum(ps_g_analytical)  # normalize

avg_cases_gen5[3] = average_cases(gen5_target_20per)

avg_cases_gen10[3] = average_cases(gen10_target_20per)



####
# Random at 15%

psi_g = np.loadtxt(f"nerccs2022/geo_05-5-22_random_rollout_15_5_intv.txt", delimiter=',')
inverted_s_m = psi_g.T
ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
random_constant_gen5 = ps_g_analytical / np.sum(ps_g_analytical)  # normalize

psi_g = np.loadtxt(f"nerccs2022/geo_05-5-22_random_rollout_15_10_intv.txt", delimiter=',')
inverted_s_m = psi_g.T
ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
random_constant_gen10 = ps_g_analytical / np.sum(ps_g_analytical)  # normalize

rand_gen5_line = average_cases(random_constant_gen5)

rand_gen10_line = average_cases(random_constant_gen10)



axs[0,0].set_ylabel("Number of cases")
axs[0,0].set_title("Average cases at Gen 5")
axs[0,1].set_title("Average cases at Gen 10")

axs[0,0].scatter(case_vals, avg_cases_gen5, color = color[0])
axs[0,0].axhline(rand_gen5_line, linestyle = 'dashed', color = color[0])

axs[0,1].scatter(case_vals, avg_cases_gen10, color = color[0])
axs[0,1].axhline(rand_gen10_line, linestyle = 'dashed', color = color[0])

# PLOTTING THE WORST CASE SCENARIO 

worst_gen5 = np.zeros((points))

worst_gen10 = np.zeros((points))

thresh = 10**(-3)

rand_gen5_line = worst_case_scenario(random_constant_gen5, thresh)

rand_gen10_line = worst_case_scenario(random_constant_gen10, thresh)


worst_gen5[0] = worst_case_scenario(gen5_target_2point5per, thresh)

worst_gen10[0] = worst_case_scenario(gen10_target_2point5per, thresh)


worst_gen5[1] = worst_case_scenario(gen5_target_5per, thresh)

worst_gen10[1] = worst_case_scenario(gen10_target_5per, thresh)

worst_gen5[2] = worst_case_scenario(gen5_target_10per, thresh)

worst_gen10[2] = worst_case_scenario(gen10_target_10per, thresh)

worst_gen5[3] = worst_case_scenario(gen5_target_20per, thresh)

worst_gen10[3] = worst_case_scenario(gen10_target_20per, thresh)



axs[1,0].set_ylabel("Number of cases")
axs[1,0].set_title("Worst case at Gen 5")
axs[1,1].set_title("Worst case at Gen 10")


axs[1,0].scatter(case_vals, worst_gen5, color = color[1])
axs[1,0].axhline(rand_gen5_line, linestyle = 'dashed', color = color[1])

axs[1,1].scatter(case_vals, worst_gen10, color = color[1])
axs[1,1].axhline(rand_gen10_line, linestyle = 'dashed', color = color[1])


# PLOTTING SHIT HITS THE FAN

fan_gen5 = np.zeros((points))

fan_gen10 = np.zeros((points))

thresh = 300

rand_gen5_line = s_hits_fan(random_constant_gen5, thresh)

rand_gen10_line = s_hits_fan(random_constant_gen10, thresh)


fan_gen5[0] = s_hits_fan(gen5_target_2point5per, thresh)

fan_gen10[0] = s_hits_fan(gen10_target_2point5per, thresh)


fan_gen5[1] = s_hits_fan(gen5_target_5per, thresh)

fan_gen10[1] = s_hits_fan(gen10_target_5per, thresh)

fan_gen5[2] = s_hits_fan(gen5_target_10per, thresh)

fan_gen10[2] = s_hits_fan(gen10_target_10per, thresh)

fan_gen5[3] = s_hits_fan(gen5_target_20per, thresh)

fan_gen10[3] =s_hits_fan(gen10_target_20per, thresh)




axs[2,0].set_ylabel("Probability")
axs[2,0].semilogy()
axs[2,0].set_title("S**t hits the fan at Gen 5")
axs[2,1].set_title("S**t hits the fan at Gen 10")



axs[2,0].scatter(case_vals, fan_gen5, color = color[2])
axs[2,0].axhline(rand_gen5_line, linestyle = 'dashed', color = color[2])

axs[2,1].scatter(case_vals, fan_gen10, color = color[2])
axs[2,1].axhline(rand_gen10_line, linestyle = 'dashed', color = color[2])

# PLOTTING THE POINTLESS INTERVENTION

pointless_gen5 = np.zeros((points))

pointless_gen10 = np.zeros((points))

thresh = 399


psi_g = np.loadtxt(f"nerccs2022/geo_05-5-22_non_interv_5.txt", delimiter=',') 

inverted_s_m = psi_g.T
ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
gen5_none = ps_g_analytical / np.sum(ps_g_analytical)  # normalize

psi_g = np.loadtxt(f"nerccs2022/geo_05-5-22_non_interv_10.txt", delimiter=',') 

inverted_s_m = psi_g.T
ps_g_analytical = np.sum(inverted_s_m, axis=0) #marginalize over s
gen10_none = ps_g_analytical / np.sum(ps_g_analytical)  # normalize


rand_gen5_line = pointless_interv(random_constant_gen5, gen5_none, thresh)



rand_gen10_line = pointless_interv(random_constant_gen10, gen10_none, thresh)


pointless_gen5[0] = pointless_interv(gen5_target_2point5per, gen5_none, thresh)

pointless_gen10[0] = pointless_interv(gen10_target_2point5per, gen10_none, thresh)


pointless_gen5[1] = pointless_interv(gen5_target_5per, gen5_none, thresh)

pointless_gen10[1] = pointless_interv(gen10_target_5per, gen10_none, thresh)

pointless_gen5[2] = pointless_interv(gen5_target_10per,  gen5_none, thresh)

pointless_gen10[2] = pointless_interv(gen10_target_10per, gen10_none, thresh)

pointless_gen5[3] = pointless_interv(gen5_target_20per,  gen5_none, thresh)

pointless_gen10[3] = pointless_interv(gen10_target_20per, gen10_none, thresh)



axs[3,0].set_ylabel("Probability")
axs[3,0].set_title("Pointless intervention at Gen 5")
axs[3,1].set_title("Pointless intervention at Gen 10")

axs[3,0].semilogy()


axs[3,0].set_xlim(0,25)
axs[3,1].set_xlim(0,25)
axs[3,0].set_xlabel("% Vaccinated for each rollout")
axs[3,1].set_xlabel("% Vaccinated for each rollout")

axs[3,0].scatter(case_vals, pointless_gen5, color = color[3])
axs[3,0].axhline(rand_gen5_line, linestyle = 'dashed', color = color[3])

axs[3,1].scatter(case_vals, pointless_gen10, color = color[3])
axs[3,1].axhline(rand_gen10_line, linestyle = 'dashed', color = color[3])
