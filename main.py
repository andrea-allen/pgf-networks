import matplotlib.pyplot as plt
import networkx as nx
import random
import numpy as np
import itertools
import scipy
from scipy import stats
from matplotlib import rc
import math


def infections_caused_matrix(P_k, beta, x=1):
    dist = np.zeros((len(P_k), len(P_k)))
    for k in range(len(P_k)):
        for l in range(k):
            dist[l][k] = P_k[k]*p_l_infected(k, l, beta)*(x**l)
    return dist

def p_l_infected(k, l, beta):
    return scipy.special.binom(k, l)*(beta**l)*((1-beta)**(k-l))

def little_test():
    degree_data = np.random.poisson(10, 100)
    P_k = pdf_of(degree_data)
    plt.plot(P_k)
    plt.show()
    g_0 = infections_caused_matrix(P_k, .15)
    avg_degree = z1_of(P_k)
    g_1 = g1_of(P_k)
    return g_0

def pdf_of(degree_list):
    g0 = np.zeros(len(degree_list))
    for i in range(len(g0)):
        g0[i] = len(np.where(degree_list == i)[0]) / len(degree_list)
    return g0

def z1_of(g_0):
    z1=0
    for k in range(len(g_0)):
        z1+=(k * g_0[k])
    return z1

def g1_of(g_0):
    g_1 = np.zeros(len(g_0))
    for k in range(len(g_0)-1):
        g_1[k] = (k+1)*g_0[k+1]
    return g_1/(z1_of(g_0))



# h = 100;
#
# for lambda =1: h
#
# maxk = h;
#
# x0 = 0.5;
#
# g1 = poisspdf(1:maxk, lambda ); % Poisson distribution set up
#
#
#
#                              g1 = g1./ sum(g1); % normalizes the data for G1

#                              g0 =[0 g1]; % gives the polynomial constants in G0
#
#                              % which is the PGF for the prob distribution
#
#                              % of vertex degrees k.
#
# for k=2: maxk
#
# g0(k) = g1(k - 1) / (k - 1);
#
# end
#
# g0 = g0. / sum(g0);
#
# fun =
#
#
# @(x)
#
#
# polyval(fliplr(g1), x) - x;
#
# u = fzero(fun, x0);
#
# S(lambda ) = 1-polyval(fliplr(g0), u);
#
# end
#
# plot(1: 10, S(1: 10))


if __name__ == '__main__':
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    rc('text', usetex=True)
    print('pgfs yay!')
    little_test()