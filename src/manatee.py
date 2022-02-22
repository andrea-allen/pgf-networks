#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
On defining the most absurd tansmission expression ever.
"""

import unittest
import math


## Ancillary functions 


# This function yields the index of the most recent intervention to k (w function)
# Updated as of 1/31 - not using this function
def index_of_recent_intervene(k, interGen):
    for i in range(1, len(interGen)):
        # print(i)
        if (k >= interGen[i - 1]) & (k < interGen[i]):
            if i - 1 == 0:
                return 0
            else:
                return i - 1

    return len(interGen) - 1


# This function gives the next soonest intervention (f function)
# Updated as of 1/31 - not using this function
def next_soonest_intervene(k, interGen):
    for i in range(1, len(interGen)):
        if (k >= interGen[i - 1]) & (k < interGen[i]):
            return interGen[i]

    return math.inf  # infinity


# Takes terms out of the equation when k is not an intervention gen 
# Updated as of 1/31 - not using this function
def indicator_1(k):
    if k == math.inf:
        return 0

    return 1


# Designates when k is one of the intervention gens
# def indicator_2(k, interGen):

#     if k in interGen:
#         return 1
#     return 0


# Amount of time from k until generation next_soonest_intervene
# Updated as of 1/31 - not using this function
def time_from_current_to_next_intervene(k, l, interGen):
    f_l = next_soonest_intervene(l, interGen)

    if f_l - k > 0:
        return f_l - k
    return -1


# Poisson Process Time term
# Updated as of 1/31 - not using this function
def pp_time_term(b, q, gen, f_gen, interGen):
    return q * b / (time_from_current_to_next_intervene(gen, f_gen, interGen))


# Function for defining parameters
# def define_parameters(beta_vec, gamma_vec, q_vec, vaccs_vec, gen, interGen):
#     beta = beta_vec[index_of_recent_intervene(gen, interGen)]
#     gamma = gamma_vec[index_of_recent_intervene(gen, interGen)]
#     q = q_vec[index_of_recent_intervene(gen, interGen)]
#     vacc = vaccs_vec[index_of_recent_intervene(gen, interGen)]
#
#     return [beta, gamma, q, vacc]


def first_term_transmission(betas, gammas, qs, vaccs, gen, interGen):
    params = define_parameters(betas, gammas, qs, vaccs, gen, interGen)

    beta_first = params[0]
    gamma_first = params[1]
    q_first = params[2]
    vacc_first = params[3]

    numerator_first = beta_first * (1 - vacc_first)
    denom_first = beta_first + gamma_first + indicator_1(next_soonest_intervene(gen, interGen)) * (
            1 - vacc_first) * pp_time_term(beta_first, q_first, gen, gen, interGen)

    return numerator_first / denom_first


def transmission_sum_term(betas, gammas, qs, vaccs, gen, interGen):
    params_g = define_parameters(betas, gammas, qs, vaccs, gen, interGen)

    beta_sum_g = params_g[0]
    gamma_sum_g = params_g[1]
    q_sum_g = params_g[2]
    vacc_sum_g = params_g[3]

    first_sum_term = (1 - vacc_sum_g) * pp_time_term(beta_sum_g, q_sum_g, gen, gen, interGen) \
                     / (beta_sum_g + gamma_sum_g
                        + indicator_1(next_soonest_intervene(gen, interGen))
                        * (1 - vacc_sum_g) * pp_time_term(beta_sum_g, q_sum_g, gen, gen, interGen))

    working_sum = 0

    for ell in range(gen + 1, max(interGen) + 1):

        if ell in interGen:
            params_l = define_parameters(betas, gammas, qs, vaccs, ell, interGen)

            beta_sum_l = params_l[0]
            gamma_sum_l = params_l[1]
            q_sum_l = params_l[2]
            vacc_sum_l = params_l[3]

            second_sum_term = beta_sum_l * (1 - vacc_sum_l) / (beta_sum_l + gamma_sum_l
                                                               + indicator_1(next_soonest_intervene(ell, interGen))
                                                               * (1 - vacc_sum_l) * pp_time_term(beta_sum_l, q_sum_l,
                                                                                                 gen, ell, interGen))

            working_sum += first_sum_term * second_sum_term

    return working_sum


def transmission_expression(betas, gammas, qs, vaccs, gen, interGen):
    first_term = first_term_transmission(betas, gammas, qs, vaccs, gen, interGen)
    second_term = transmission_sum_term(betas, gammas, qs, vaccs, gen, interGen)

    return first_term + second_term


def l_of_g(g, l, b_lminus, b_l, gamma_lminus, gamma_l, q_lminus, q_l, v_lminus, v_l):
    first_term = ((1 - v_lminus) * q_lminus * b_lminus) / (
                b_lminus + gamma_lminus + b_lminus * q_lminus * (1 - v_lminus))
    second_term = ((1 - v_l) * b_l) / (b_l + gamma_l + b_l * q_l * (1 - v_l))

    return first_term ** (l - g) * second_term


def t_of_g(betas, gammas, qs, vaccs, g):
    g_term = ((1 - vaccs[g]) * betas[g]) / (betas[g] + gammas[g] + betas[g] * qs[g] * (1 - vaccs[g]))
    # .266
    betas_extended = list(betas)
    vaccs_extended = list(vaccs)
    gammas_extended = list(gammas)
    qs_extended = list(qs)
    for i in range(500):
        betas_extended.append(betas[-1])
        vaccs_extended.append(vaccs[-1])
        gammas_extended.append(gammas[-1])
        qs_extended.append(qs[-1])
    for l in range(g + 1, len(betas_extended)):
        # g_term += l_of_g(g, l, betas[l - 1], betas[l], gammas[l - 1], gammas[l], qs[l - 1], qs[l], vaccs[l - 1],
        #                  vaccs[l])
        g_term += l_of_g(g, l, betas_extended[l - 1], betas_extended[l], gammas_extended[l - 1], gammas_extended[l],
                         qs_extended[l - 1], qs_extended[l], vaccs_extended[l - 1],
                         vaccs_extended[l])

    return g_term
