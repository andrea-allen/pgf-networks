import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime

# Say we use a serial interval of 5 days
# then g is given by... 0, 5, 10, 15, 20 days intervals after first case
# plot data point: for each state, (s) cumulative cases at g=20 (serial days after first recorded case)
# write a function that's: get first case date
# another function to get the date 20 days after first case (iloc + 20)
# then can plot each state for g=20 and (s)

#    Next, vary g and start on a first date. Then for each g along the x axis, get the corresponding multiple of the
# serial interval. Plot the number of cases at that serial interval after the first case.

def run():
    print('Covid-19')
    colors = ['red', 'orange', 'blue', 'yellow', 'purple']
    covid_df = read_jh_data()
    all_states = list(pd.unique(covid_df['Province_State']))
    y_labels = []
    for state in ['California', 'New York', 'Texas', 'Washington', 'Vermont']:
        #TODO convert to datetime
        state_df = get_by_state(covid_df, state)
        state_df_summed = sum_for_state(state_df)
        first_case = get_first_case_date(state_df_summed)
        start_date = '3/15/20'
        start_date_as_date = datetime.datetime.strptime(start_date, '%m/%d/%y')
        for g in range(1, 50):
            cases_20_days = add_serial_interval(start_date, 5*g)
            if cases_20_days[0] =='0':
                cases_20_days = cases_20_days[1:]
            if cases_20_days[2] == '0':
                cases_20_days = cases_20_days[0:2] + cases_20_days[3:]
            if cases_20_days[2] == '/' and cases_20_days[3] == '0':
                cases_20_days = cases_20_days[0:3] + cases_20_days[4:]
            # s_cum_first = state_df_summed.loc[first_case]
            s_cum = state_df_summed.loc[cases_20_days]
            cases_20_days = datetime.datetime.strptime(cases_20_days, '%m/%d/%y').date()
            # first_case = datetime.datetime.strptime(first_case, '%m/%d/%y').date()
            plt.scatter(g, s_cum, label=state, color=colors[0])
            if g%5==0:
                plt.text(g, s_cum, f'{state}, {cases_20_days} \n {s_cum} cases')
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
    return covid_df

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

