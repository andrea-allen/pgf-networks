#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
On defining the most absurd tansmission expression ever.
"""

import unittest
import math


## Ancillary functions 


# This function yields the index of the most recent intervention to k (w function)
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
def next_soonest_intervene(k, interGen):
    for i in range(1, len(interGen)):
        if (k >= interGen[i - 1]) & (k < interGen[i]):
            return interGen[i]

    return math.inf  # infinity


# Takes terms out of the equation when k is not an intervention gen 
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
def time_from_current_to_next_intervene(k, l, interGen):
    f_l = next_soonest_intervene(l, interGen)

    if f_l - k > 0:
        return f_l - k
    return -1


# Poisson Process Time term
def pp_time_term(b, q, gen, f_gen, interGen):
    return q * b / (time_from_current_to_next_intervene(gen, f_gen, interGen))


# Function for defining parameters
def define_parameters(beta_vec, gamma_vec, q_vec, vaccs_vec, gen, interGen):
    beta = beta_vec[index_of_recent_intervene(gen, interGen)]
    gamma = gamma_vec[index_of_recent_intervene(gen, interGen)]
    q = q_vec[index_of_recent_intervene(gen, interGen)]
    vacc = vaccs_vec[index_of_recent_intervene(gen, interGen)]

    return [beta, gamma, q, vacc]


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


# %% TEST


class TestTransmissionFuncts(unittest.TestCase):

    def test_index_of_recent_intervene(self):
        h = [0, 3, 5]
        self.assertEqual(index_of_recent_intervene(1, h), 0)
        self.assertEqual(index_of_recent_intervene(3, h), 1)
        self.assertEqual(index_of_recent_intervene(4, h), 1)

    def test_next_soonest_intervene(self):
        h = [0, 3, 5]
        self.assertEqual(next_soonest_intervene(1, h), 3)
        self.assertEqual(next_soonest_intervene(3, h), 5)
        self.assertEqual(next_soonest_intervene(4, h), 5)
        self.assertEqual(next_soonest_intervene(5, h), math.inf)

    def test_time_from_current_to_next_intervene(self):
        h = [0, 3, 5]
        self.assertEqual(time_from_current_to_next_intervene(1, 1, h), 2)
        self.assertEqual(time_from_current_to_next_intervene(3, 3, h), 2)
        self.assertEqual(time_from_current_to_next_intervene(4, 4, h), 1)
        self.assertEqual(time_from_current_to_next_intervene(5, 5, h), math.inf)

        self.assertEqual(time_from_current_to_next_intervene(1, 2, h), 2)
        self.assertEqual(time_from_current_to_next_intervene(3, 4, h), 2)
        self.assertEqual(time_from_current_to_next_intervene(4, 6, h), math.inf)
        self.assertEqual(time_from_current_to_next_intervene(5, 6, h), math.inf)

    def test_first_term_transmission_basecase(self):
        h = [0]
        beta = [0.9]
        gamma = [.001]
        expected_T = 0.9 / (0.9 + 0.001)
        q = [2.5]
        vacc = [0]
        self.assertEqual(first_term_transmission(beta, gamma, q, vacc, 0, h), expected_T)
        self.assertEqual(first_term_transmission(beta, gamma, q, vacc, 1, h), expected_T)
        self.assertEqual(first_term_transmission(beta, gamma, q, vacc, 3, h), expected_T)

    def test_first_term_transmission_trivial(self):
        h = [0]
        beta = [0.9]
        gamma = [.001]
        expected_T = 0.9 / (0.9 + 0.001)
        q = [2.5]
        vacc = [0]
        self.assertEqual(first_term_transmission(beta, gamma, q, vacc, 0, h), expected_T)
        self.assertEqual(first_term_transmission(beta, gamma, q, vacc, 1, h), expected_T)
        self.assertEqual(first_term_transmission(beta, gamma, q, vacc, 3, h), expected_T)

    def test_transmission_expression_basecase(self):
        h = [0]
        beta = [0.9]
        gamma = [.001]
        expected_T = 0.9 / (0.9 + 0.001)
        q = [2.5]
        vacc = [0]
        self.assertEqual(transmission_expression(beta, gamma, q, vacc, 0, h), expected_T)
        self.assertEqual(transmission_expression(beta, gamma, q, vacc, 1, h), expected_T)
        self.assertEqual(transmission_expression(beta, gamma, q, vacc, 2, h), expected_T)
        self.assertEqual(transmission_expression(beta, gamma, q, vacc, 3, h), expected_T)

    def test_transmission_expression_trivial(self):
        h = [0, 2]
        beta = [0.9, 0.9]
        gamma = [.001, .001]
        expected_T = 0.9 / (0.9 + 0.001)
        q = [2.5, 2.5]
        vacc = [0, 0]
        self.assertEqual(transmission_expression(beta, gamma, q, vacc, 0, h), expected_T)
        self.assertEqual(transmission_expression(beta, gamma, q, vacc, 1, h), expected_T)
        self.assertEqual(transmission_expression(beta, gamma, q, vacc, 2, h), expected_T)
        self.assertEqual(transmission_expression(beta, gamma, q, vacc, 3, h), expected_T)

    def test_transmission_expression_vaccination_change_single(self):
        h = [0, 2]
        beta = [0.9, 0.9]
        gamma = [.001, .001]
        expected_T0 = (.9) / (.901 + (2.5 * .9) / (2)) + (
                ((2.5 * .9) / (2)) / (.901 + (2.5 * .9) / (2)) * ((.9 * .5) / (.901)))
        expected_T1 = (.9) / (.901 + (2.5 * .9) / (1)) + (
                ((2.5 * .9) / (1)) / (.901 + (2.5 * .9) / (1)) * ((.9 * .5) / (.901)))
        expected_T2 = (.9 * .5) / (.901)
        q = [2.5, 2.5]
        vacc = [0, .5]
        self.assertEqual(transmission_expression(beta, gamma, q, vacc, 0, h), expected_T0)
        self.assertEqual(transmission_expression(beta, gamma, q, vacc, 1, h), expected_T1)
        self.assertEqual(transmission_expression(beta, gamma, q, vacc, 2, h), expected_T2)
        self.assertEqual(transmission_expression(beta, gamma, q, vacc, 3, h), expected_T2)

    def test_transmission_expression_vaccination_change_double(self):
        h = [0, 2, 4]
        beta = [0.9, 0.9, 0.9]
        gamma = [0.001, 0.001, 0.001]

        expected_T0 = ((.9) / (.901 + ((1 * 2.5 * 0.9) / 2))) \
                      + ((((1 * 2.5 * 0.9) / 2) / (.901 + (1 * 2.5 * 0.9) / 2))
                         * (0.9 * .5) / (.901 + (.5 * 2.5 * 0.9) / 4)) \
                      + ((((2.5 * .9) / 2) / (.901 + (2.5 * .9) / 2))
                         * ((.9 * .3) / .901))

        expected_T1 = ((.9) / (.901 + ((2.5 * 0.9) / 1))) \
                      + ((((1 * 2.5 * 0.9) / 1) / (.901 + (1 * 2.5 * 0.9) / 1))
                         * (0.9 * .5) / (.901 + (.5 * 2.5 * 0.9) / 3)) \
                      + ((((1 * 2.5 * .9) / 1) / (.901 + (1 * 2.5 * .9) / 1))
                         * ((.9 * .3) / .901))

        expected_T2 = ((.9 * .5) / (.901 + ((.5 * 2.5 * 0.9) / 2))) \
                      + ((((.5 * 2.5 * 0.9) / 2) / (.901 + (.5 * 2.5 * 0.9) / 2))
                         * ((0.9 * .3) / (.901)))

        expected_T3 = ((.9 * .5) / (.901 + ((.5 * 2.5 * 0.9) / 1))) \
                      + ((((.5 * 2.5 * 0.9) / 1) / (.901 + (.5 * 2.5 * 0.9) / 1))
                         * ((0.9 * .3) / (.901)))

        expected_T4 = (.9 * .3) / (.901)
        q = [2.5, 2.5, 2.5]
        vacc = [0, .5, .7]  # cumulative numbers
        self.assertEqual(transmission_expression(beta, gamma, q, vacc, 0, h), expected_T0)  # good
        self.assertAlmostEqual(transmission_expression(beta, gamma, q, vacc, 1, h), expected_T1, 5)  # good
        self.assertAlmostEqual(transmission_expression(beta, gamma, q, vacc, 2, h), expected_T2, 5)  # good
        self.assertEqual(transmission_expression(beta, gamma, q, vacc, 3, h), expected_T3)  # good
        self.assertAlmostEqual(transmission_expression(beta, gamma, q, vacc, 4, h), expected_T4, 5)  # good

    def test_transmission_expression_triple(self):
        h = [0, 3, 4, 6]
        beta = [.8, .8, .8, .8]
        beta_const = .8
        gamma_const = 0.2
        gamma = [.2, .2, .2, .2]
        q = [2.5, 2.5, 2.5, 2.5]
        q_const = 2.5
        vacc_cum = [0, .2, .4, .7] ## cumulative random vax

        expected_T0 = None #TODO
        expected_T1 = (beta_const*(1-0))\
                      /(beta_const+gamma_const + (q_const*beta_const/2)) \
                      + ( ((q_const*beta_const/2)/(beta_const+gamma_const + (q_const*beta_const/2)))
           *(beta_const*(1-.2)/(beta_const+gamma_const+((1-.2)*q_const*beta_const/3))) ) \
                      + (((q_const * beta_const / 2) / (beta_const + gamma_const + (q_const * beta_const / 2)))
           * (beta_const * (1 - .4) / (beta_const + gamma_const + ((1 - .4) * q_const * beta_const / 5))))\
                      + (((q_const * beta_const / 2) / (beta_const + gamma_const + (q_const * beta_const / 2)))
           * (beta_const * (1 - .7) / (beta_const + gamma_const )))

        expected_T2 = None #TODO
        expected_T3 = None #TODO
        expected_T4 = None #TODO
        print(f'\n{expected_T1}')
        self.assertAlmostEqual(transmission_expression(beta, gamma, q, vacc_cum, 1, h), expected_T1, 5)  # good





if __name__ == '__main__':
    unittest.main()
