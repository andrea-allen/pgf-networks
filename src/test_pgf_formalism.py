#%% TEST FORMALISM

import unittest
from src.pgf_formalism import *
import numpy as np


class TestFormalismChanges(unittest.TestCase):

    def test_H_target(self): # Good
        degree_d = np.array([0.4, 0.3, 0.2, 0.1])
        delta_k = np.zeros((1, len(degree_d)))
        prop = 0.5
        g = 0

        k_crit, temp_prop, delta_g_k = critical_degree_calc(prop, degree_d, delta_k, g)

        H_g_code = np.sum([(k+1) * delta_g_k[g][k+1] * degree_d[k+1] for k in range(0, len(degree_d)-1)]) /np.sum([(k+1) * degree_d[k+1] for k in range(0, len(degree_d)-1)])

        print(H_g_code)
        true_H = ( (1*delta_g_k[g][1]*degree_d[1]) + (2 * delta_g_k[g][2] * degree_d[2]) + (3 * delta_g_k[g][3] * degree_d[3]) )/( (1 * degree_d[1]) + (2 * degree_d[2]) + (3 * degree_d[3]) )
        print(true_H)
        self.assertAlmostEqual(H_g_code, true_H,5)  # good


    def test_gofg(self):
        degree_d = np.array([0.4, 0.3, 0.2, 0.1])
        delta_k = np.zeros((1, len(degree_d)))
        prop = 0.5
        g = 0

        k_crit, temp_prop, delta_g_k = critical_degree_calc(prop, degree_d, delta_k, g)

        updated = modify_g0(degree_d, k_crit, temp_prop)
        expected_g_g = np.array([0,0,0])
        code_g_g = g_g_of(updated, delta_g_k, g)
        self.assertEqual(code_g_g.all(), expected_g_g.all())

    def test_g0_of(self):
        degree_d = np.array([0.4, 0.3, 0.2, 0.1])
        delta_k = np.zeros((1, len(degree_d)))
        prop = 0.5
        g = 0

        k_crit, temp_prop, delta_g_k = critical_degree_calc(prop, degree_d, delta_k, g)

        new_g = modify_g0(degree_d, k_crit, temp_prop)

        #print(new_g)

        correct_change = np.array([0.4, 0.1, 0, 0])
        correct_change = correct_change/np.sum(correct_change)
        #print(correct_change)

        self.assertEqual(new_g.all(), correct_change.all())

    def test_critical_degree_deltas(self): #Good
        degree_d = np.array([0.4, 0.3, 0.2, 0.1])
        delta_k = np.zeros((1, len(degree_d)))
        prop = 0.5
        g = 0

        k_crit, temp_prop, expected_array = critical_degree_calc(prop, degree_d, delta_k, g)
        print(k_crit)
        print(temp_prop)

        final_array = np.array([0, degree_d[k_crit]*0.1, 1, 1])
        print(final_array)
        print(expected_array)
        self.assertEqual(expected_array.all(), final_array.all())



if __name__ == '__main__':
    unittest.main()

