import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import pandas as pd
import datetime
from src import pgf_formalism
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
import matplotlib.contour as ctr
import seaborn as sns

# Say we use a serial interval of 5 days
# then g is given by... 0, 5, 10, 15, 20 days intervals after first case
# plot data point: for each state, (s) cumulative cases at g=20 (serial days after first recorded case)
# write a function that's: get first case date
# another function to get the date 20 days after first case (iloc + 20)
# then can plot each state for g=20 and (s)

#    Next, vary g and start on a first date. Then for each g along the x axis, get the corresponding multiple of the
# serial interval. Plot the number of cases at that serial interval after the first case.

def exploratory_state_plot(state_name, show_state=False):
    covid_df = read_jh_data()
    state_df = get_by_state(covid_df, state_name)
    state_df_summed = sum_for_state(state_df)
    first_case = get_first_case_date(state_df_summed)
    cases_20_days = add_serial_interval(first_case, 30)
    if cases_20_days[0] == '0':
        cases_20_days = cases_20_days[1:]
    if cases_20_days[2] == '0':
        cases_20_days = cases_20_days[0:2] + cases_20_days[3:]
    if cases_20_days[2] == '/' and cases_20_days[3] == '0':
        cases_20_days = cases_20_days[0:3] + cases_20_days[4:]

    early_cases = pd.DataFrame(state_df_summed.loc[first_case:cases_20_days])
    early_cases['dates'] = pd.to_datetime(early_cases.index)
    plt.plot(early_cases['dates'], early_cases[0], label=state_name)
    plt.xticks(rotation=45)
    if show_state:
        plt.show()

def covid_data(ax, cp, fontcolor, solo_plot=True):
    print('Covid-19')
    colors = ['red', 'orange', 'blue', 'yellow', 'purple']
    covid_df = read_jh_data()
    # all_states = list(pd.unique(covid_df['Province_State']))
    y_labels = []
    state_labels = {'California': 'CA','Washington':'WA','Vermont':'VT','New Mexico':'NM','Michigan':'MI','Wyoming':'WY',
                    'New York':'NY', 'Arizona':'AZ', 'Illinois':'IL', 'Louisiana':'LA', 'South Dakota':'SD', 'Hawaii':'HI',
                    'Massachusetts':'MA', 'Georgia':'GA'}
    state_number = 0
    for state in ['California','Washington','Illinois','Hawaii', 'Massachusetts', 'Georgia', 'Arizona']:
        state_number += 1
        #TODO convert to datetime
        state_df = get_by_state(covid_df, state)
        state_df_summed = sum_for_state(state_df)
        first_case = get_first_case_date(state_df_summed)
        # start_date = '3/15/20'
        start_date = first_case
        # start_date = state_df_summed.head(1).index.values[0]
        start_date_as_date = datetime.datetime.strptime(start_date, '%m/%d/%y')
        line_plot_vals_s_cum = []
        line_plot_vals_g = []
        for g in range(1, 12):
            cases_20_days = add_serial_interval(start_date, 4*g) #starts at g=1, where start_date is g=0
            if cases_20_days[0] == '0':
                cases_20_days = cases_20_days[1:]
            if cases_20_days[2] == '0':
                cases_20_days = cases_20_days[0:2] + cases_20_days[3:]
            if cases_20_days[2] == '/' and cases_20_days[3] == '0':
                cases_20_days = cases_20_days[0:3] + cases_20_days[4:]
            # s_cum_first = state_df_summed.loc[first_case]
            if g == 1:
                gen_1_date = cases_20_days
            last_index_date = cases_20_days
            s_cum = state_df_summed.loc[cases_20_days]
            cases_20_days = datetime.datetime.strptime(cases_20_days, '%m/%d/%y').date()
            # first_case = datetime.datetime.strptime(first_case, '%m/%d/%y').date()
            # ax.annotate(f'{state_labels[state]}, {cases_20_days} \n {s_cum} cases', (cx, cy), color='black', weight='bold', fontsize=10, ha='center', va='center')
            # ax.scatter(g-1 + .5, s_cum, label=state_labels[state], color='navy', s=.4) # Plotting offset by 1 since graph starts at g=1, not g=0
            line_plot_vals_g.append(g-1)
            line_plot_vals_s_cum.append(s_cum)
            if g % 1 == 0:
                if s_cum > 10 and g <= 6:
                    this_color = 'black' #change if necessary
                elif s_cum > 80:
                    this_color = 'black'
                else:
                    this_color = 'black'
                # Plotting offset by 1 (g-1 instead of g) since graph starts at g=1, not g=0
                if s_cum==1:
                    cases = 'case'
                else:
                    cases = 'cases'
                if g==1:
                    # Plotting first dates in upper left, since they all get squished below
                    plt.text(g-1, 1000-10*state_number, f'{state_labels[state]}\n {cases_20_days}', fontsize=10, color=this_color, fontweight='bold') # Plotting offset by 1 since graph starts at g=1, not g=0
                elif g==11:
                    state_date_final_dict = {'HI':'Apr 20', 'MA': 'Mar 13', 'CA': 'Mar 10', 'WA':'Mar 6', 'IL':'Mar 8', 'AZ':'Mar 10', 'GA':'offmap'}
                    plt.text(g - 1, s_cum, f'{state_labels[state]}\n {state_date_final_dict[state_labels[state]]}', fontsize=10,
                             color=this_color, fontweight='bold')  # Plotting offset by 1 since graph starts at g=1, not g=0
                else:
                    plt.text(g - 1, s_cum, f'{state_labels[state]}',
                             fontsize=10, color=this_color, fontweight='bold')  # Plotting offset by 1 since graph starts at g=1, not g=0
                # ctr.ClabelText(g, s_cum, f'{state_labels[state]}, {cases_20_days} \n {s_cum} cases')
            # plt.scatter(first_case, s_cum_first, label=state, color=colors[0])
            # plt.text(first_case + datetime.timedelta(days=1), s_cum_first, f'{state} \n {first_case} \n {s_cum_first} cases')
            # y_labels.append(s_cum_first)
            y_labels.append(s_cum)
        state_data = state_df_summed.loc[first_case:last_index_date]
        # state_data_shifted = np.zeros(len(state_data))
        # for s in range(len(state_data)):
        #     case_count = state_data[s]
        #     if case_count < 10:
        #         state_data_shifted[s] = case_count + 0.5
        #     elif case_count < 100:
        #         state_data_shifted[s] = case_count + 5
        #     elif case_count < 1000:
        #         state_data_shifted[s] = case_count + 50
        #     elif case_count < 2000:
        #         state_data_shifted[s] = case_count + 500
        #     else:
        #         state_data_shifted[s] = case_count + 1000
        ax.plot(np.arange(0, 11, 11/len(np.array(state_data))), np.array(state_data), color='black', lw=1)
        ax.plot(np.arange(0, 11, 11/len(np.array(state_data))), np.array(state_data), color='white', lw=6, alpha=0.3)
        # colors.remove(colors[0])
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

    # fmt = {}
    # strs = ['first', 'second', 'third', 'fourth', 'fifth', 'sixth', 'seventh']
    # for l, s in zip(cp.levels[10:17], strs):
    #     fmt[l] = s

    # Label every other level using strings
    # ax.clabel(cp, cp.levels[10:17], inline=True, fmt=fmt, fontsize=10)

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
    g = 9
    s = 18

    k_vals = np.arange(.09, 1.1, .2)
    r0_vals = np.arange(1.2, 5, .2)
    Z_vals = np.zeros((len(k_vals), len(r0_vals)))
    for i in range(len(r0_vals)):
        for j in range(len(k_vals)):
            r0 = r0_vals[i]
            k = k_vals[j]
            g0, g1 = pgf_formalism.offspring_dists(r0=r0, k=k, p0=0.03, length=500) #TODO pull derivation of p0 from contour file from LHD
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

    # np.savetxt('X_vals_g9_s18.txt', X)
    # np.savetxt('Y_vals_g9_s18.txt', Y)
    # np.savetxt('Z_vals_g9_s18.txt', Z_vals)

    return X, Y, Z_vals

"""
Helper function used for making Figure 3 contour plot
"""
def make_figure2():
    k = .1
    r0 = 2.5
    g_vals = np.arange(1, 13)
    s_vals = np.arange(0, 1001)
    Z_vals = np.zeros((len(s_vals), len(g_vals)))
    p0 = 1 - (k * ((k / (k + r0)) ** (k - 1) - 1) / (1 - k))
    g0, g1 = pgf_formalism.offspring_dists(r0=r0, k=k, p0=0, length=1001)  # TODO pull derivation of p0 from contour file from LHD
    results = pgf_formalism.compute_extinct_prob_all(n_gens=13, renorm=True, custom_g0=g0, custom_g1=g1)
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
    for i in range(len(g_vals)):
        for s in range(len(s_vals)):
            if s <= i+1:
                Z_vals[s][i] = 0

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

def contour_figs(data_path=None):
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



    if data_path is None:
        X, Y, Z = make_figure()
    else:
        X = np.loadtxt(f'{data_path}/X_vals2.txt') ## PRESAVED VALUES
        Y = np.loadtxt(f'{data_path}/Y_vals2.txt')
        Z = np.loadtxt(f'{data_path}/Z_vals2.txt')
    pal = sns.color_palette("colorblind")
    print(pal.as_hex())
    cp = ax.contour(X, Y, Z,levels=40, cmap=get_custom_cmap())
    # cbar = fig.colorbar(cp)
    ax.clabel(cp, inline=True, fontsize=10)
    # Patch:
    # ax.add_patch(plt.Rectangle((1.4,0.1), 2.5, 0.54, fill=False,
    #                            edgecolor="r", linewidth=3))

    ax.set_ylabel(r'Dispersion parameter $k$')
    ax.set_xlabel(r'$R_0$')
    ax.set_xticks(np.arange(1.2, 5, .5))
    ax.set_yticks([1/10, 2/10, 3/10, 4/10, 5/10, 6/10, 7/10, 8/10, 9/10, 10/10])
    ax.set_yticklabels(['.1', '.2', '.3', '.4', '.5', '.6', '.7', '.8', '.9', '1.0'])

    # cbar.set_label(r'Probability of epidemic survival', rotation=270)

    norm = mpl.colors.Normalize(vmin=cp.cvalues.min(), vmax=cp.cvalues.max())
    # a previous version of this used
    # norm= matplotlib.colors.Normalize(vmin=cs.vmin, vmax=cs.vmax)
    # which does not work any more
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cp.cmap)
    sm.set_array([])
    fig.colorbar(sm, ticks=np.arange(min(cp.levels), max(cp.levels), .1), label=r'Probability of epidemic survival')
    # cbar.set_label(r'Probability of epidemic survival', rotation=270)

    # Save to file.
    plt.tight_layout(0.1)
    # plt.savefig('contour_draft2.png')
    # plt.savefig('ccontour_WA_9_18.svg', format='svg')
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
    cp = ax.contour(X, Y, Z,levels=40)
    covid_data(ax, cp, solo_plot=False)
    # cbar = fig.colorbar(cp)
    # ax.clabel(cp, inline=True, fontsize=10)
    # Patch:
    # ax.add_patch(plt.Rectangle((1.4,0.1), 2.5, 0.54, fill=False,
    #                            edgecolor="r", linewidth=3))

    ax.set_ylabel(r'Cumulative cases $s$')
    ax.set_xlabel(r'Epidemic generation $g$')
    # ax.set_xlim([0, 11])
    # cbar.set_label(r'Probability of epidemic survival', rotation=270)
    norm = mpl.colors.Normalize(vmin=cp.cvalues.min(), vmax=cp.cvalues.max())
    # a previous version of this used
    # norm= matplotlib.colors.Normalize(vmin=cs.vmin, vmax=cs.vmax)
    # which does not work any more
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cp.cmap)
    sm.set_array([])
    fig.colorbar(sm, ticks=np.arange(min(cp.levels), max(cp.levels), .1), label=r'Probability of epidemic survival')

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
    # plt.savefig('contour2_svg2.svg', format='svg')

    plt.show()

"""
Figure 3 here is what is used in Paper 1
"""
def contour_fig3():
    plt.rcParams["text.usetex"] = True
    plt.rcParams["font.size"] = 16
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Fira Sans", "PT Sans", "Open Sans", "Roboto", "DejaVu Sans", "Liberation Sans",
                                       "sans-serif"]
    plt.rcParams["xtick.major.width"] = 2
    plt.rcParams["xtick.major.size"] = 8
    plt.rcParams["ytick.major.width"] = 2
    plt.rcParams["ytick.major.size"] = 8

    plt.rcParams.update({'font.size': 12})

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
    ax = plt.gca()
    # fig, (ax) = plt.subplots(1, 1, figsize=(16, 6))
    # fig, ((ax)) = plt.subplots(1,1)
    # fig.subplots_adjust(bottom=0.15)

    ax.set_yscale('log')
    # ax.set_ylim(5,100)

    X, Y, Z = make_figure2()
    # TODO don't plot s=0 row
    Z = Z[1:]

    # cp = ax.imshow(Z, interpolation='none', extent=extent, origin='lower', alpha=1, aspect='auto',
    #                cmap=get_custom_cmap()) #Spectral_r
    cp = ax.imshow(Z, interpolation='none', origin='lower', alpha=1, aspect='auto',
                   cmap=get_custom_cmap()) #Spectral_r
    # ax.yaxis.grid(True, which='minor')
    # ax.xaxis.grid(True, which='minor')
    ax.set_xticks(np.arange(0, 12))
    ax.set_xticklabels(np.arange(1,13))
    ytick_spots = 10**(np.array([0, np.log10(5), 1, np.log10(50), 2, np.log10(500), 3]))
    ax.set_yticks(ytick_spots)
    ax.set_yticklabels(np.array([1, 5, 10, 50, 100, 500, 1000]))
    ax.tick_params(axis='y', which='minor', left=False)
    # cp = ax.contour(X, Y, Z,levels=40)
    covid_data(ax, cp, fontcolor='white',solo_plot=False)
    ax.set_ylim(0.75, 1000)
    ax.set_xlim(0, 11.5)
    cbar = fig.colorbar(cp)
    # ax.clabel(cp, inline=True, fontsize=10)
    # Patch:
    # ax.add_patch(plt.Rectangle((1.4,0.1), 2.5, 0.54, fill=False,
    #                            edgecolor="r", linewidth=3))

    ax.set_ylabel(r'Cumulative cases $s$', fontsize=14)
    ax.set_xlabel(r'Epidemic generation $g$', fontsize=14)
    # ax.set_xlim([0, 11])
    cbar.set_label(r'Probability of epidemic survival', rotation=270, fontsize=14)
    # norm = mpl.colors.Normalize(vmin=cp.cvalues.min(), vmax=cp.cvalues.max())
    # a previous version of this used
    # norm= matplotlib.colors.Normalize(vmin=cs.vmin, vmax=cs.vmax)
    # which does not work any more
    # sm = plt.cm.ScalarMappable(norm=norm, cmap=cp.cmap)
    # sm.set_array([])
    # fig.colorbar(sm, ticks=np.arange(min(cp.levels), max(cp.levels), .1), label=r'Probability of epidemic survival')

    # plt.text(0.05, 0.05, r'low $R_0$, high variance', fontsize=16, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    # plt.text(0.57, 0.95, r'high $R_0$, low variance', fontsize=16, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)


    # Save to file.
    # plt.tight_layout(0.1)
    # plt.savefig('contour3_yscale_fixed1.svg', format='svg')

    plt.show()

def get_custom_cmap():
    hex_list = ['#0173b2', '#56b4e9',  '#029e73', '#ece133', '#de8f05', '#d55e00', '#cc78bc', '#fbafe4', '#ca9161',  '#949494',  ]
    hex_list_less = ['#0173b2', '#56b4e9',  '#ece133', '#de8f05', '#d55e00',] #yellow: removed
    # my_cmap = get_continuous_cmap(hex_list)
    # x, y = np.mgrid[-5:5:0.05, -5:5:0.05]
    # z = (np.sqrt(x ** 2 + y ** 2) + np.sin(x ** 2 + y ** 2))
    # # hex_list = ['#0091ad', '#3fcdda', '#83f9f8', '#d6f6eb', '#fdf1d2', '#f8eaad', '#faaaae', '#ff57bb']
    #
    # fig, ax = plt.subplots(1, 1)
    # im = ax.imshow(z, cmap=get_continuous_cmap(hex_list_less))
    # fig.colorbar(im)
    # ax.yaxis.set_major_locator(plt.NullLocator())  # remove y axis ticks
    # ax.xaxis.set_major_locator(plt.NullLocator())  # remove x axis ticks
    # plt.show()
    return get_continuous_cmap(hex_list_less)


def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.

        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0, 1, len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]






