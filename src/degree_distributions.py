import numpy as np
import math


def power_law_degree_distrb(maxk=40):
    p_k = np.empty(maxk)
    p_k[0] = 0
    for k in range(1, maxk):
        p_k[k] = (k ** (-2)) * (math.e ** (-k / 5))
    p_k = p_k / np.sum(p_k)
    return p_k


def binomial_degree_distb(N, lam=6):
    p_k = np.empty(40)
    p = lam / N
    for k in range(0, len(p_k)):
        p_k[k] = (p ** k) * ((1 - p) ** (N - k)) * math.comb(N, k)
    return p_k
