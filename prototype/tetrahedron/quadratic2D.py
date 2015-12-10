import numpy as np

# Matrix to obtain q_i from e-values on the input simplex
# (q_i) = M(e_i)
# q_0 = epsilon_0,0
# q_1 = epsilon_1,0
# q_2 = epsilon_0,1
# q_3 = epsilon_1/2,0
# q_4 = epsilon_1/2,1/2
# q_5 = epsilon_0,1/2
G = [[  1,  0,  0,  0,  0,  0],
     [ -3, -1,  0,  4,  0,  0],
     [ -3,  0, -1,  0,  0,  4],
     [  2,  2,  0, -4,  0,  0],
     [  4,  0,  0, -4,  4, -4],
     [  2,  0,  2,  0,  0, -4]]

def get_mesh2D(mesh):
    return np.array([[x[1], x[0]] for x in np.ndindex((mesh[0], mesh[1]))])

class QuadraticTetrahedron2D:
    def __init__(self, e, prec=1e-5):
        self._q = np.dot(G, e)
        self._prec = prec

    def show_mesh(self, mesh):
        print(get_mesh2D(mesh))

    def run(self):
        self._to_ellipsoid()

    def _to_ellipsoid(self):
        A1 = self._to_ellipsoid_step1()
        A2 = self._to_ellipsoid_step2()
        t3 = self._to_ellipsoid_step3()

    def _to_ellipsoid_step1(self):
        q = self._q
        S = 1.0 / 2 * np.array([[2 * q[3], q[4]], [q[4], 2 * q[5]]])
        w, A = np.linalg.eigh(S) # S = U^T.diag(w).U
        t = q[[1, 2]]
        q[1:3] = np.dot(A.T, t)
        q[4] = 0
        q[3], q[5] = w
        return A

    def _to_ellipsoid_step2(self):
        q = self._q
        if (abs(q[3]) < self._prec and
            abs(q[5]) > self._prec):
            q[3] = q[5]
            q[5] = 0
            A = np.array([[0, 1], [1, 0]], dtype=float)
            t = q[[1, 2]]
            q[1:3] = np.dot(A.T, t)
            return A
        else:
            return np.eye(2)
        
    def _to_ellipsoid_step3(self):
        q = self._q
        if (abs(q[3]) > self._prec and
            abs(q[1]) > self._prec):
            return 
        else:
            return np.zeros(2)
        
