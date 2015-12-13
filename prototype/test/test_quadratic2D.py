import unittest
import sys
import numpy as np
from tetrahedron.quadratic2D import QuadraticTetrahedron2D

def _one(x, y):
    return 1

def _cos(x, y, n=1):
    return np.cos(n * np.pi * x) + cos(n * np.pi * y)

class TestQuadraticTetrahedron2D(unittest.TestCase):
    def setUp(self):
        e = np.arange(6, dtype='double')
        mesh = [4, 4]
        self._qt2d = QuadraticTetrahedron2D(e, mesh)

    def tearDown(self):
        pass
    
    def test_run(self):
        self._qt2d.run()

    def test_show_mesh(self):
        self._qt2d.show_mesh()

if __name__ == '__main__':
    unittest.main()
