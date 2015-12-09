import numpy as np

# Matrix to obtain q_i from e-values on the input simplex
# (q_i) = M(e_i)
# q_0 = epsilon_0,0
# q_1 = epsilon_1,0
# q_2 = epsilon_0,1
# q_3 = epsilon_1/2,0
# q_4 = epsilon_1/2,1/2
# q_5 = epsilon_0,1/2
M = [[  1,  0,  0,  0,  0,  0],
     [ -3, -1,  0,  4,  0,  0],
     [ -3,  0, -1,  0,  0,  4],
     [  2,  2,  0, -4,  0,  0],
     [  4,  0,  0, -4,  4, -4],
     [  2,  0,  2,  0,  0, -4]]


def get_mesh2D(mesh):
    return np.array([[x[1], x[0]] for x in np.ndindex((mesh[0], mesh[1]))])

class QuadraticTetrahedron2D:
    def __init__(self, e):
        self._q = np.dot(M, e)

    def show_mesh(self, mesh):
        print(get_mesh2D(mesh))

    def show_M(self):
        print(np.array(M))

    def run(self):
        self._to_ellipsoid()

    def _to_ellipsoid(self):
        q = self._q
        r = q[0]
        t = q[1:3]
        S = 1.0 / 2 * np.array([[2 * q[3], q[4]], [q[4], 2 * q[5]]])
        w, U = np.linalg.eigh(S)
        print(S)
        print(w)
        print(U)
        print(np.dot(np.dot(U.T, np.diag(w)), U))
