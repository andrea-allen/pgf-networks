#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy.stats import binom
"""
Created on Tue Apr 13 15:11:09 2021

@author: nicholasjroberts
"""

# takes a deg dist and transmission prob and returns offspring dist
def gen_offspring_dist(deg_dist, T, custom=False):
    if not custom:
        d_dist = np.zeros(len(deg_dist))
        d_dist[0:-1] = np.multiply(deg_dist[1:], range(1,len(deg_dist)))
        d_dist = d_dist / deg_dist.dot(range(len(deg_dist)))
    else:
        d_dist = deg_dist
    l_dist = np.zeros(len(d_dist))
    for l in range(len(l_dist)):
        for k in range(l,len(l_dist)):
            l_dist[l] += binom.pmf(l,k,T)*d_dist[k]
    return l_dist/np.sum(l_dist)

def G1(ddist):
    new_dist = np.zeros(len(ddist))
    new_dist[0:-1] = np.multiply(ddist[1:],range(1,len(ddist)))
    new_dist = new_dist / ddist.dot(range(len(ddist)))
    
    Pk = np.zeros((len(new_dist), 2))
    Pk[:,0] = range(len(new_dist))
    Pk[:,1] = list(new_dist)
    return Pk

def G(x, pk):
    return np.power(x[np.newaxis].T, pk[:,0]).dot(pk[:,1])
    #return np.multiply(np.power(x, pk[:,0]), pk[:,1])

def GT(T, pk):
    N = pk.shape[0]+20
    z = np.exp(2 * np.pi * complex(0,1) * np.arange(N) / N)
    G_at_x = G(z*T+1-T, pk)
    return np.roll(np.absolute(np.fft.ifft(G_at_x))[::-1],1)

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
    
def sngl_ext_prob(deg_dist, T, fft=True, custom=False):
    #return ext_prob_iter(make_pgf(gen_offspring_dist(deg_dist, T)))
    if fft and custom:
        return ext_prob_iter( make_pgf( GT(T, deg_dist) ) )
    elif fft:
        return ext_prob_iter( make_pgf( GT(T, G1(deg_dist)) ) )
    else:
        return ext_prob_iter(make_pgf(gen_offspring_dist(deg_dist, T, custom)))

"""
for Psi cake indexed as [g][s][m]
compute extinction probability array
"""
def gen_ext_prob_array(psi, d_dist, T, fft=True, custom=False):
    e_prob = sngl_ext_prob(d_dist, T, fft, custom)
    extnct_array = np.zeros(psi.shape)
    for g in range(extnct_array.shape[0]):
        for s in range(extnct_array.shape[1]):
            for m in range(extnct_array.shape[2]):
                extnct_array[g][s][m] = psi[g][s][m]*e_prob**m
    return extnct_array