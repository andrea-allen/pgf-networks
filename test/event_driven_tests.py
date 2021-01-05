import unittest


class TestNode(unittest.TestCase):
    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

    def test_failing_test(self):
        self.assertEqual('2', '4')


class TestEdge(unittest.TestCase):
    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

    def test_failing_test(self):
        self.assertEqual('2', '4')


if __name__ == '__main__':
    unittest.main()
