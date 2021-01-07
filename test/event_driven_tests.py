import unittest
import numpy as np
from src.event_driven import *
from src.SIR_sims import *

class TestGeneral(unittest.TestCase):
    def test_array_fills(self):
        my_zero_array = np.zeros((10, 10))
        for i in range(10):
            for j in range(10):
                my_zero_array[i][j] = 3.1
        my_full_array = np.full((10, 10), 3.1)
        self.assertTrue((my_full_array == my_zero_array).all())

class TestNode(unittest.TestCase):
    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')


class TestEdge(unittest.TestCase):
    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

class TestSimulation(unittest.TestCase):
    def setUp(self):
        N = 100
        self.beta = 0.99
        self.gamma = 0.01
        dd = power_law_degree_distrb()
        graph, pos = generate_graph(N, dd)
        self.Lambda = np.full((N, N), self.beta)
        self.Gamma = np.full(N, self.gamma)
        self.graph = graph

    def test_IS_edge_is_removed_after_single_step(self):
        simulation = Simulation(100000, self.graph, self.beta, self.gamma, None)
        simulation.intialize()
        starting_IS_list_length = len(simulation.V_IS)
        self.assertGreaterEqual(starting_IS_list_length, 1)
        self.assertEqual(len(simulation.V_I), 1)
        simulation.single_step()
        self.assertEqual(len(simulation.V_I), 2)
        self.assertEqual(len(simulation.V_IS), starting_IS_list_length-1)
        current_length_IS_list = len(simulation.V_IS)
        if current_length_IS_list > 0:
            simulation.single_step()
            self.assertEqual(len(simulation.V_I), 3)
            self.assertEqual(len(simulation.V_IS), current_length_IS_list - 1)

    def test_intervention(self):
        simulation = Simulation(100000, self.graph, self.beta, self.gamma, None)
        simulation.intialize()
        starting_IS_list_length = len(simulation.V_IS)
        self.assertGreaterEqual(starting_IS_list_length, 1)
        self.assertEqual(len(simulation.V_I), 1)
        self.assertEqual(simulation.V_IS[0].event_rate, self.beta)
        new_beta = 0.8
        simulation.intervene(new_beta)
        self.assertGreaterEqual(starting_IS_list_length, 1)
        self.assertEqual(len(simulation.V_I), 1)
        self.assertEqual(simulation.V_IS[0].event_rate, new_beta)


if __name__ == '__main__':
    unittest.main()
