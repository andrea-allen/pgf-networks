#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy.stats import binom
"""
Created on Tue Apr 13 15:11:09 2021

@author: nicholasjroberts
"""

# takes a deg dist and transmission prob and returns offspring dist
def gen_offspring_dist(deg_dist, T):
    l_dist = np.zeros(len(deg_dist))
    for l in range(len(l_dist)):
        for k in range(l,len(l_dist)):
            l_dist[l] += binom.pmf(l,k,T)*deg_dist[k]
    return l_dist/np.sum(l_dist)

# takes a distribution and returns a pgf as a function
def make_pgf(dist):
    def ret_pgf(x):
        ret = 0
        for k in range(len(dist)):
            ret += dist[k]*x**k
        return ret
    return ret_pgf 

"""
a recursive method for finding extinction probability
mathematically it is function composition
call should go as ext_prob_iter(offspring_pgf)
the arg mu is the offspring pgf
"""
def ext_prob_iter(mu, x=0, comp=1, tol=10**(-7)):
    if np.abs(mu(x) - comp) < tol:
        return mu(x)
    else:
        return ext_prob_iter(mu, mu(x), mu(x), tol)
    
def sngl_ext_prob(deg_dist, T):
    return ext_prob_iter(make_pgf(gen_offspring_dist(deg_dist, T)))

"""
for Psi cake indexed as [g][s][m]
compute extinction probability array
"""
def gen_ext_prob_array(psi, d_dist, T):
    e_prob = sngl_ext_prob(d_dist, T)
    extnct_array = np.zeros(psi.shape)
    for g in range(extnct_array.shape[0]):
        for s in range(extnct_array.shape[1]):
            for m in range(extnct_array.shape[2]):
                extnct_array[g][s][m] = psi[g][s][m]*e_prob**m
    return extnct_array