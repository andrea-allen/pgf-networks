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
    degree_data = np.random.poisson(10, 50)
    P_k = pdf_of(degree_data)
    plt.plot(P_k)
    plt.show()
    g_0_matrix = infections_caused_matrix(P_k, .15)
    avg_degree = z1_of(P_k)
    z2 = z1_of(g1_of(g1_of(P_k)))   # This might be wrong
    beta = .3
    r_0 = beta*z2/avg_degree
    g_1 = g1_of(P_k)
    g_1_matrix = infections_caused_matrix(g_1, .15)
    # myplot = plt.plot(g_0_matrix)
    fig, ax = plt.subplots()
    fig1 = ax.imshow(g_0_matrix, cmap='plasma')
    plt.show()
    fig, ax = plt.subplots()
    fig2 = ax.imshow(g_1_matrix, cmap='plasma')
    plt.show()
    Psi_sm = phase_space(P_k, g_1)
    fig, ax = plt.subplots()
    fig3 = ax.imshow(Psi_sm[0], cmap='plasma')
    plt.show()
    return g_0_matrix

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

def z2_of(g_0):
    z2=0
    for k in range(len(g_0)-2):
        z2+=(k+1)*(k+2)*g_0[k+2]
    return z2


def g1_of(g_0):
    g_1 = np.zeros(len(g_0))
    for k in range(len(g_0)-1):
        g_1[k] = (k+1)*g_0[k+1]
    return g_1/(z1_of(g_0))

def phase_space(g_0, g_1, g=10):
    Psi_sm = np.zeros((10, 100,100))
    # Initial condition:
    Psi_sm[0][1][1] = 1
    return Psi_sm

# Given a degree distribution for G0 (or the degree distribution of the entire network.

# degreeDist = insert probability list for the starting distribution
# maxk = maximum k
# start_G0 = pdf_of(degreeDist) Maybe adjust this to have the generation in this and
    # loop for as many generations as we see fit.
#
# start_G1 = g1_of(g0)
#

# for g in range(len(generation)):
#     begin taking away nodes from the dist that become infected
#       Then recompute G1 for each time step until hitting intervention
#       Different function comes in there where we decide on the new distribution or new T
#       For the generation before intervention, solve the selfconsistent equation and
#       and compute final outbreak size. 
#

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