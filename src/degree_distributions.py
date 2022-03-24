import numpy as np
import math
import src.pgf_formalism as pgf


def power_law_degree_distrb(maxk=40, alpha=2, mu=5):
    p_k = np.empty(maxk)
    p_k[0] = 0
    for k in range(1, maxk):
        p_k[k] = (k ** (-alpha)) * (math.e ** (-k / mu))
    p_k = p_k / np.sum(p_k)
    return p_k


def binomial_degree_distb(N, lam=6):
    p_k = np.empty(N)
    p = lam / N
    for k in range(0, len(p_k)):
        try:
            p_k[k] = (p ** k) * ((1 - p) ** (N - k)) * math.comb(N, k)
        except OverflowError:
            p_k[k] = 0 #why zero? Because for large k, p**k will go to zero
            print(N, k)
    return p_k

def exponential_degree_dist(N, p):
    p_k = np.empty(N)
    for k in range(1, N):
        p_k[k] =  (1 - p) * p**k
    p_k = p_k/np.sum(p_k)
    return p_k

def mean_degree(deg_distr):
    total_sum = 0
    for k in range(len(deg_distr)):
        total_sum += k*deg_distr[k]
    return total_sum

def chain_degree_dist(N):
    p_k = np.zeros(N)
    p_k[0] = 0
    p_k[1] = .01
    p_k[2] = .99
    # p_k[3] = .29
    # p_k[4] = .3
    return p_k




def compute_T_threshold_powerlaw(deg_distr, mu):
    # Currently only works for power law networks
    k_expected = mean_degree(deg_distr)
    k_square_expected = math.e**(2/mu) / (math.e**(1/mu) - 1)
    thresh = k_expected / (k_square_expected - k_expected)
    return thresh

def compute_T_threshold_binom(deg_distr):
    # Currently only works for binomial networks
    k_expected = mean_degree(deg_distr)
    pgf_double_prime = np.zeros(len(deg_distr))
    g1 = np.zeros(len(deg_distr))
    for k in range(len(deg_distr) - 1):
        g1[k] = (k + 1) * deg_distr[k + 1]
    for j in range(len(g1) - 1):
        pgf_double_prime[j] = (j+1)*g1[j+1]
    k_square_expected = np.sum(pgf_double_prime) + k_expected # var_k + mean_k ^2
    thresh = k_expected / (k_square_expected - k_expected)
    return thresh
