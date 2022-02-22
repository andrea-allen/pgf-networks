# %% TEST
import unittest
from src.manatee import *
import numpy as np

class TestTransmissionFuncts(unittest.TestCase):

    # def test_index_of_recent_intervene(self):
    #     h = [0, 3, 5]
    #     self.assertEqual(index_of_recent_intervene(1, h), 0)
    #     self.assertEqual(index_of_recent_intervene(3, h), 1)
    #     self.assertEqual(index_of_recent_intervene(4, h), 1)
    #
    # def test_next_soonest_intervene(self):
    #     h = [0, 3, 5]
    #     self.assertEqual(next_soonest_intervene(1, h), 3)
    #     self.assertEqual(next_soonest_intervene(3, h), 5)
    #     self.assertEqual(next_soonest_intervene(4, h), 5)
    #     self.assertEqual(next_soonest_intervene(5, h), math.inf)
    #
    # def test_time_from_current_to_next_intervene(self):
    #     h = [0, 3, 5]
    #     self.assertEqual(time_from_current_to_next_intervene(1, 1, h), 2)
    #     self.assertEqual(time_from_current_to_next_intervene(3, 3, h), 2)
    #     self.assertEqual(time_from_current_to_next_intervene(4, 4, h), 1)
    #     self.assertEqual(time_from_current_to_next_intervene(5, 5, h), math.inf)
    #
    #     self.assertEqual(time_from_current_to_next_intervene(1, 2, h), 2)
    #     self.assertEqual(time_from_current_to_next_intervene(3, 4, h), 2)
    #     self.assertEqual(time_from_current_to_next_intervene(4, 6, h), math.inf)
    #     self.assertEqual(time_from_current_to_next_intervene(5, 6, h), math.inf)
    #
    # def test_first_term_transmission_basecase(self):
    #     h = [0]
    #     beta = [0.9]
    #     gamma = [.001]
    #     expected_T = 0.9 / (0.9 + 0.001)
    #     q = [2.5]
    #     vacc = [0]
    #     self.assertEqual(first_term_transmission(beta, gamma, q, vacc, 0, h), expected_T)
    #     self.assertEqual(first_term_transmission(beta, gamma, q, vacc, 1, h), expected_T)
    #     self.assertEqual(first_term_transmission(beta, gamma, q, vacc, 3, h), expected_T)

    # def test_first_term_transmission_trivial(self):
    #     h = [0]
    #     beta = [0.9]
    #     gamma = [.001]
    #     expected_T = 0.9 / (0.9 + 0.001)
    #     q = [2.5]
    #     vacc = [0]
    #     self.assertEqual(first_term_transmission(beta, gamma, q, vacc, 0, h), expected_T)
    #     self.assertEqual(first_term_transmission(beta, gamma, q, vacc, 1, h), expected_T)
    #     self.assertEqual(first_term_transmission(beta, gamma, q, vacc, 3, h), expected_T)

    def test_transmission_expression_basecase(self):
        beta = [0.9, 0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9]
        gamma = [0.001, 0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001, 0.001]
        q = [2.5, 2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5]
        vacc = [0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        expected_T = 0.9 / (0.9 + 0.001)
        # THIS DOESNT STAY SAME
        # self.assertTrue(t_of_g(beta, gamma, q, vacc, 0)==t_of_g(beta, gamma, q, vacc, 1))
        self.assertEqual(t_of_g(beta, gamma, q, vacc, 0), t_of_g(beta, gamma, q, vacc, 1))
        # self.assertTrue(t_of_g(beta, gamma, q, vacc, 1)==t_of_g(beta, gamma, q, vacc, 2))
        # self.assertTrue(t_of_g(beta, gamma, q, vacc, 2)==t_of_g(beta, gamma, q, vacc, 3))
        # self.assertEqual(expected_T, t_of_g(beta, gamma, q, vacc, 0))
        # self.assertEqual(expected_T, t_of_g(beta, gamma, q, vacc, 1))
        # self.assertEqual(expected_T, t_of_g(beta, gamma, q, vacc, 2))
        # self.assertEqual(expected_T, t_of_g(beta, gamma, q, vacc, 3))

    def test_transmission_expression_vaccination_change_single(self):
        # single intervention of 50% random at generation 2
        len_gens = 20
        beta = np.full(len_gens, .9)
        gamma = np.full(len_gens, .001)
        q = np.full(len_gens, 2.5)
        vacc = list(np.full(len_gens, 0))
        vacc[3:] = list(np.full(len_gens-3, 0.5)) # At gen 2, 50% random vaccination enacted
        # TODO T0 and T1, T2 should get changed
        expected_T0 = (.9) / (.901 + (2.5 * .9) / (2)) + (
                ((2.5 * .9) / (2)) / (.901 + (2.5 * .9) / (2)) * ((.9 * .5) / (.901)))
        expected_T1 = (.9) / (.901 + (2.5 * .9) / (1)) + (
                ((2.5 * .9) / (1)) / (.901 + (2.5 * .9) / (1)) * ((.9 * .5) / (.901)))
        expected_T2 = (.9 * .5) / (.901) #TODO change
        expected_T3 = (.9 * .5) / (.901)
        expected_T4 = (.9 * .5) / (.901)
        print(vacc)
        # self.assertEqual(expected_T0, t_of_g(beta, gamma, q, vacc, 0))
        # self.assertEqual(expected_T1, t_of_g(beta, gamma, q, vacc, 1))
        self.assertEqual(expected_T3, t_of_g(beta, gamma, q, vacc, 3))
        self.assertEqual(expected_T4, t_of_g(beta, gamma, q, vacc, 4))

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