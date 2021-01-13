import unittest
import numpy as np
from src.event_driven import *
from src.SIR_sims import *
import networkx as nx
import timeit
import time

class TestGeneral(unittest.TestCase):
    def test_array_fills(self):
        my_zero_array = np.zeros((10, 10))
        for i in range(10):
            for j in range(10):
                my_zero_array[i][j] = 3.1
        my_full_array = np.full((10, 10), 3.1)
        self.assertTrue((my_full_array == my_zero_array).all())

class TestNode(unittest.TestCase):
    def test_node_attributes(self):
        node = Node(24, 0, 0, 0.01)
        self.assertEqual(node.label, 24)
        self.assertEqual(node.gen, 0)
        self.assertEqual(node.state, 0)
        self.assertEqual(node.event_rate, 0.01)

        node.infect()
        self.assertEqual(node.state, 1)

        node.recover()
        self.assertEqual(node.state, 2)

        node.set_gen(3)
        self.assertEqual(node.gen, 3)

        other_node = Node(32, 0, 0, 0.01)
        same_node = Node(24, 0, 0, 0.01)
        same_node_different_state = Node(24, 1, 1, 0.01)
        self.assertTrue(node.equals(same_node))
        self.assertTrue(node.equals(same_node_different_state))
        self.assertFalse(node.equals(other_node))


class TestEdge(unittest.TestCase):
    def test_edge_attributes(self):
        infected_node = Node(24, 0, 1, 0.01)
        susceptible_node = Node(32, -1, 0, 0.01)
        IS_edge = Edge(infected_node, susceptible_node, 0.8)
        self.assertEqual(IS_edge.i.label, 24)
        self.assertEqual(IS_edge.j.label, 32)
        self.assertEqual(IS_edge.j.state, 0)
        IS_edge.infect()
        self.assertEqual(IS_edge.j.state, 1)
        self.assertEqual(IS_edge.j.gen, 1)
        self.assertEqual(IS_edge.i.state, 1)
        self.assertEqual(IS_edge.i.gen, 0)

class TestSimulation(unittest.TestCase):
    def setUp(self):
        N = 10
        self.beta = 0.99
        self.gamma = 0.0001
        dd = power_law_degree_distrb()
        graph, pos = generate_graph(N, dd)
        graph = nx.generators.complete_graph(N)
        N = len(graph.nodes())
        self.Lambda = np.full((N, N), self.beta)
        self.Gamma = np.full(N, self.gamma)
        self.graph = graph

    def test_IS_edges_are_updated_after_single_step(self):
        adjacency_matrix = np.array(nx.adjacency_matrix(self.graph).todense())
        simulation = Simulation(100000, self.graph, self.beta, self.gamma, self.Lambda, self.Gamma, None,
                                adjacency_matrix)
        simulation.intialize()
        print('before single step')
        simulation.display_info()
        self.assertGreaterEqual(len(simulation.V_IS), 1)
        self.assertEqual(len(simulation.V_I), 1)
        simulation.single_step()
        print('after single step')
        self.assertEqual(len(simulation.V_I), 2)
        simulation.display_info()
        infected_nodes = list(map(lambda node: node.label, simulation.V_I))
        for edge in simulation.V_IS:
            self.assertIn(edge.i.label, infected_nodes)
            self.assertNotIn(edge.j.label, infected_nodes)
        if len(simulation.V_IS) > 0:
            simulation.single_step()
            print('after second single step')
            simulation.display_info()
            self.assertEqual(len(simulation.V_I), 3)
            infected_nodes = list(map(lambda node: node.label, simulation.V_I))
            for edge in simulation.V_IS:
                self.assertIn(edge.i.label, infected_nodes)
                self.assertNotIn(edge.j.label, infected_nodes)


    def test_intervention(self):
        adjacency_matrix = np.array(nx.adjacency_matrix(self.graph).todense())
        simulation = Simulation(100000, self.graph, self.beta, self.gamma, self.Lambda, self.Gamma, None, adjacency_matrix)
        simulation.intialize()
        starting_IS_list_length = len(simulation.V_IS)
        self.assertGreaterEqual(starting_IS_list_length, 1)
        self.assertEqual(len(simulation.V_I), 1)
        self.assertEqual(simulation.V_IS[0].event_rate, self.beta)
        new_beta = 0.8
        simulation.intervene(new_beta, True)
        self.assertGreaterEqual(starting_IS_list_length, 1)
        self.assertEqual(len(simulation.V_I), 1)
        self.assertEqual(simulation.V_IS[0].event_rate, new_beta)

        simulation = Simulation(100000, self.graph, self.beta, self.gamma, self.Lambda, self.Gamma, None,
                                adjacency_matrix)
        simulation.intialize()
        starting_IS_list_length = len(simulation.V_IS)
        self.assertGreaterEqual(starting_IS_list_length, 1)
        self.assertEqual(len(simulation.V_I), 1)
        self.assertEqual(simulation.V_IS[0].event_rate, self.beta)
        new_beta = 0.8
        simulation.intervene(new_beta, False)
        self.assertGreaterEqual(starting_IS_list_length, 1)
        self.assertEqual(len(simulation.V_I), 1)
        self.assertEqual(simulation.V_IS[0].event_rate, self.beta)

    #TODO write a test that tests if intervention happens at the correct time, need to modify the code structure for this

    @unittest.skip('Work in progress')
    def test_intervention_occurs_correctly(self):
        #TODO use mocks
        # assert that when intervene is called, only one generation member of that gen is exists
        self.assertFalse(True)


if __name__ == '__main__':
    unittest.main()
