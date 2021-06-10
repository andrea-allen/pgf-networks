import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import scipy.stats as spt
import src.gen_extinct_prob
from src import pgf_formalism
from matplotlib.colors import LinearSegmentedColormap

# Say we use a serial interval of 5 days
# then g is given by... 0, 5, 10, 15, 20 days intervals after first case
# plot data point: for each state, (s) cumulative cases at g=20 (serial days after first recorded case)
# write a function that's: get first case date
# another function to get the date 20 days after first case (iloc + 20)
# then can plot each state for g=20 and (s)

#    Next, vary g and start on a first date. Then for each g along the x axis, get the corresponding multiple of the
# serial interval. Plot the number of cases at that serial interval after the first case.

def covid_data(ax, solo_plot=True):
    print('Covid-19')
    colors = ['red', 'orange', 'blue', 'yellow', 'purple']
    covid_df = read_jh_data()
    # all_states = list(pd.unique(covid_df['Province_State']))
    y_labels = []
    for state in ['Michigan', 'New York', 'Wyoming', 'California']:
        #TODO convert to datetime
        state_df = get_by_state(covid_df, state)
        state_df_summed = sum_for_state(state_df)
        first_case = get_first_case_date(state_df_summed)
        start_date = '3/15/20'
        start_date = first_case
        # start_date = state_df_summed.head(1).index.values[0]
        start_date_as_date = datetime.datetime.strptime(start_date, '%m/%d/%y')
        for g in range(1, 5):
            cases_20_days = add_serial_interval(start_date, 3*g)
            if cases_20_days[0] == '0':
                cases_20_days = cases_20_days[1:]
            if cases_20_days[2] == '0':
                cases_20_days = cases_20_days[0:2] + cases_20_days[3:]
            if cases_20_days[2] == '/' and cases_20_days[3] == '0':
                cases_20_days = cases_20_days[0:3] + cases_20_days[4:]
            # s_cum_first = state_df_summed.loc[first_case]
            s_cum = state_df_summed.loc[cases_20_days]
            cases_20_days = datetime.datetime.strptime(cases_20_days, '%m/%d/%y').date()
            # first_case = datetime.datetime.strptime(first_case, '%m/%d/%y').date()
            ax.scatter(g, s_cum, label=state, color=colors[0], s=.2)
            if g % 1 == 0:
                plt.text(g, s_cum, f'{state}, {cases_20_days} \n {s_cum} cases', fontsize=6)
            # plt.scatter(first_case, s_cum_first, label=state, color=colors[0])
            # plt.text(first_case + datetime.timedelta(days=1), s_cum_first, f'{state} \n {first_case} \n {s_cum_first} cases')
            # y_labels.append(s_cum_first)
            y_labels.append(s_cum)
        colors.remove(colors[0])
        # plt.plot(state_df_summed, label=state)
    # plt.xlabel(rotation=35)
    indices_to_show = np.arange(0, len(state_df_summed), 60)
    # indices_to_show = np.arange(0, 365, 30)
    # plt.xticks(state_df_summed.index.values[indices_to_show], rotation=35)
    if solo_plot:
        plt.xticks(rotation=35)
        plt.ylabel('Cumulative cases')
        plt.xlabel('Generation, (5 day increments, with serial interval=5)')
        plt.yticks(y_labels)
        plt.xticks(np.arange(0, 50))
        # plt.legend(loc='upper left')
        plt.box(on=False)
        plt.semilogy()
        # plt.text(datetime.datetime.strptime('01/30/20', '%m/%d/%y').date(), 3000, 'Serial interval: 5 days \n Points shown are date of first recorded \n case and number of cumulative \n cases, then cumulative cases 20 days later \n corresponding to generation g=4.')
        plt.show()
    # print(covid_df.head())
    return ax

def read_jh_data():
    data = pd.read_csv('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv')
    return data

def get_by_state(data, state):
    return data[data['Province_State']==state]

def sum_for_state(data):
    raw_counts = data.iloc[:,11:]
    summed_df = raw_counts.sum(axis=0)
    return summed_df

def get_first_case_date(data):
    #will return string of first date where case count is greater than 0
    first_case_date = data[data>0].head(1).index.values[0]
    return first_case_date

def add_serial_interval(date_as_str, days_to_add):
    serial_interval_days = 5
    date_as_date = datetime.datetime.strptime(date_as_str, '%m/%d/%y').date()
    date_add_days = date_as_date + datetime.timedelta(days=days_to_add)
    return date_add_days.strftime('%m/%d/%y')

def contour_fig1(r0, k, g, s):
    plt.figure('Contour plot 1')
    # Fixed g, fixed s, vary k, r0 and do the following:
    # offspring distribution, from k and r0:
    # neg_binom = spt.nbinom.pmf(n=100, p=0.5, k=2)
    neg_binom = np.arange(0, 100)
    z_vals = np.zeros((11, len(neg_binom)))
    # survival prob:
    for k in range(1, 12):
        actual_k = k/2
        # neg_binom = spt.nbinom.pmf(n=100, p=0.5, k=actual_k)
        survival_prob_matrix = np.random.rand(len(neg_binom), len(neg_binom))
        summed = np.sum(survival_prob_matrix, axis=1)
        z_vals[k-1] = summed
    x_sample = np.arange(0, 10)
    y_sample = np.arange(0, 10)
    X, Y = np.meshgrid(neg_binom, np.arange(11))
    Z1 = np.exp(-X ** 2 - Y ** 2)
    Z2 = np.exp(-(X - 1) ** 2 - (Y - 1) ** 2)
    Z = (Z1 - Z2) * 2
    # z_vals = np.full((10, 10), 2)
    # for i in range(10):
    #     z_vals[i] = np.random.randn(10)
    plt.contour(X, Y, z_vals)
    plt.show()

## Pseudocode for steps for making the contour plot
def make_figure():
    # Pick g, s as in g: 20 days into the pando, serial interval of 5, means g=4, we have s=16 cases on day 20
    # for each k we want:
    # for each r0 we want:
    # compute nbinom (g1 and g0) using k, r0
    # call the new functions we wrote to get g0, g1 with r0 and k
    # then compute entire Psi with those g1 and g0
    # compute survival probability from Psi,
    # this will be real data:
    g = 4
    s = 16

    k_vals = np.arange(.09, .25, .05)
    r0_vals = np.arange(2, 5, .2)
    Z_vals = np.zeros((len(k_vals), len(r0_vals)))
    for i in range(len(r0_vals)):
        for j in range(len(k_vals)):
            r0 = r0_vals[i]
            k = k_vals[j]
            g0, g1 = pgf_formalism.offspring_dists(r0=r0, k=k, p0=0.03, length=250) #TODO pull derivation of p0 from contour file from LHD
            results = pgf_formalism.compute_extinct_prob_all(n_gens=11, renorm=True, custom_g0=g0, custom_g1=g1)
            extnct_array = results[0] # format is g, s, m
            # then condense over m, get value at g,s
            extnct_array_all_m = extnct_array.sum(axis=2)
            val_g_s = extnct_array_all_m[g, s]
            Z_vals[j, i] = val_g_s
            # that constitutes the single value for z(r0, k)
    Z_vals = 1 - Z_vals
    X, Y = np.meshgrid(r0_vals, k_vals)
    # plt.contour(X, Y, Z_vals, levels=40)
    # plt.colorbar()
    # plt.show()

    return X, Y, Z_vals

def make_figure2():
    k = .1
    r0 = 2.5
    g_vals = np.arange(1, 10)
    s_vals = np.arange(2, 500)
    Z_vals = np.zeros((len(s_vals), len(g_vals)))
    p0 = 1 - (k * ((k / (k + r0)) ** (k - 1) - 1) / (1 - k))
    g0, g1 = pgf_formalism.offspring_dists(r0=r0, k=k, p0=0, length=1000)  # TODO pull derivation of p0 from contour file from LHD
    results = pgf_formalism.compute_extinct_prob_all(n_gens=11, renorm=True, custom_g0=g0, custom_g1=g1)
    extnct_array = results[0]  # format is g, s, m
    # then condense over m, get value at g,s
    extnct_array_all_m = extnct_array.sum(axis=2)
    for i in range(len(g_vals)):
        for j in range(len(s_vals)):
            g = g_vals[i]
            s = s_vals[j]
            val_g_s = extnct_array_all_m[g, s]
            Z_vals[j, i] = val_g_s
    Z_vals = 1 - Z_vals
    X, Y = np.meshgrid(g_vals, s_vals)
    # plt.contour(X, Y, Z_vals, levels=40)
    # plt.colorbar()
    # plt.show()

    return X, Y, Z_vals


    ## Figure 2:
    # pick r0 and k from fig 1 that are known from lit / from data / etc.
    # vary g, s
    # compute nbinom (g1 and g0) using k, r0
    # # then compute entire Psi with those g1 and g0
    # # compute survival probability from Psi,
    # # then condense over m, get value at g,s
    # populate results with each point over (g,s) survival prob
    # add real data points

def contour_figs():
    plt.rcParams["text.usetex"] = True
    plt.rcParams["font.size"] = 16
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Fira Sans", "PT Sans", "Open Sans", "Roboto", "DejaVu Sans", "Liberation Sans",
                                       "sans-serif"]
    plt.rcParams["xtick.major.width"] = 2
    plt.rcParams["xtick.major.size"] = 8
    plt.rcParams["ytick.major.width"] = 2
    plt.rcParams["ytick.major.size"] = 8

    colors = [(24 / 255, 22 / 255, 35 / 255),
              (11 / 255, 26 / 255, 69 / 255),
              (85 / 255, 114 / 255, 194 / 255),
              (216 / 255, 157 / 255, 125 / 255),
              (195 / 255, 177 / 255, 137 / 255),
              (175 / 255, 90 / 255, 59 / 255)]
    n_bins = [3, 6, 10, 100]  # Discretizes the interpolation into bins
    cmap_name = 'my_list'
    cm = LinearSegmentedColormap.from_list(
        cmap_name, colors, N=7)

    fig, ((ax)) = plt.subplots(1,1,figsize=(6,5.5), sharey=False)
    # fig.subplots_adjust(bottom=0.15)

    ax.set_yscale('log')
    # ax.set_ylim(0.005,10)

    X, Y, Z = make_figure()
    cp = ax.contour(X, Y, Z,levels=40)
    cbar = fig.colorbar(cp)
    ax.clabel(cp, inline=True, fontsize=10)
    # Patch:
    # ax.add_patch(plt.Rectangle((1.4,0.1), 2.5, 0.54, fill=False,
    #                            edgecolor="r", linewidth=3))

    ax.set_ylabel(r'Dispersion parameter $k$')
    ax.set_xlabel(r'$R_0$')
    cbar.set_label(r'Probability of epidemic survival', rotation=270)

    # plt.text(0.16, 0.33, "2019-nCoV, Wuhan",
    #          color="w",
    #          horizontalalignment='left',
    #          verticalalignment='bottom',
    #          transform=ax.transAxes,
    #          fontsize=18)
    # plt.text(0.14, 0.4, "COVID-19 (over-dispersed)",
    #          color="r",
    #          horizontalalignment='left',
    #          verticalalignment='bottom',
    #          transform=ax.transAxes,
    #          backgroundcolor="w")

    # plt.text(0.05, 0.05, r'low $R_0$, high variance', fontsize=16, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    # plt.text(0.57, 0.95, r'high $R_0$, low variance', fontsize=16, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)


    # Save to file.
    plt.tight_layout(0.1)
    plt.show()

def contour_fig2():
    plt.rcParams["text.usetex"] = True
    plt.rcParams["font.size"] = 16
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Fira Sans", "PT Sans", "Open Sans", "Roboto", "DejaVu Sans", "Liberation Sans",
                                       "sans-serif"]
    plt.rcParams["xtick.major.width"] = 2
    plt.rcParams["xtick.major.size"] = 8
    plt.rcParams["ytick.major.width"] = 2
    plt.rcParams["ytick.major.size"] = 8

    colors = [(24 / 255, 22 / 255, 35 / 255),
              (11 / 255, 26 / 255, 69 / 255),
              (85 / 255, 114 / 255, 194 / 255),
              (216 / 255, 157 / 255, 125 / 255),
              (195 / 255, 177 / 255, 137 / 255),
              (175 / 255, 90 / 255, 59 / 255)]
    n_bins = [3, 6, 10, 100]  # Discretizes the interpolation into bins
    cmap_name = 'my_list'
    cm = LinearSegmentedColormap.from_list(
        cmap_name, colors, N=7)

    fig, ((ax)) = plt.subplots(1,1,figsize=(6,5.5), sharey=False)
    # fig.subplots_adjust(bottom=0.15)

    ax.set_yscale('log')
    # ax.set_ylim(5,100)

    X, Y, Z = make_figure2()
    covid_data(ax, solo_plot=False)
    cp = ax.contour(X, Y, Z,levels=40)
    cbar = fig.colorbar(cp)
    ax.clabel(cp, inline=True, fontsize=10)
    # Patch:
    # ax.add_patch(plt.Rectangle((1.4,0.1), 2.5, 0.54, fill=False,
    #                            edgecolor="r", linewidth=3))

    ax.set_ylabel(r'Cumulative cases $s$')
    ax.set_xlabel(r'Epidemic generation $g$')
    cbar.set_label(r'Probability of epidemic survival', rotation=270)

    # plt.text(0.16, 0.33, "2019-nCoV, Wuhan",
    #          color="w",
    #          horizontalalignment='left',
    #          verticalalignment='bottom',
    #          transform=ax.transAxes,
    #          fontsize=18)
    # plt.text(0.14, 0.4, "COVID-19 (over-dispersed)",
    #          color="r",
    #          horizontalalignment='left',
    #          verticalalignment='bottom',
    #          transform=ax.transAxes,
    #          backgroundcolor="w")

    # plt.text(0.05, 0.05, r'low $R_0$, high variance', fontsize=16, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    # plt.text(0.57, 0.95, r'high $R_0$, low variance', fontsize=16, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)


    # Save to file.
    plt.tight_layout(0.1)
    plt.show()

