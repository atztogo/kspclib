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

# mesh is supported as a numpy int array.
def get_grid_address(mesh): 
    m = mesh
    return reduce_grid_address(
        m, np.array([[x[1], x[0]] for x in np.ndindex((m[0], m[1]))]))

def reduce_grid_address(mesh, address):
    return address - (address > (mesh // 2)) * mesh

def address_to_grid(mesh, address):
    a = address % mesh
    return a[1] * mesh[0] + a[0]

class QuadraticTetrahedron2D:
    def __init__(self, e, mesh, prec=1e-5):
        self._q = np.dot(G, e)
        self._mesh = np.array(mesh, dtype='intc')
        self._prec = prec

    def show_mesh(self):
        grid_address = get_grid_address(self._mesh)
        print(grid_address)
        for a in grid_address:
            print(address_to_grid(self._mesh, a))

    def run(self):
        # self._to_ellipsoid()
        pass

    def _to_ellipsoid(self):
        A1 = self._to_ellipsoid_step1()
        A2 = self._to_ellipsoid_step2()
        A3 = self._to_ellipsoid_step3()
        A4 = self._to_ellipsoid_step4()

    def _to_ellipsoid_step1(self):
        q = self._q
        S = 1.0 / 2 * np.array([[2 * q[3], q[4]], [q[4], 2 * q[5]]])
        w, v = np.linalg.eigh(S) # S = U^T.diag(w).U
        A = np.eye(3, dtype='double')
        A[0:2, 0:2] = v
        t = q[[1, 2, 0]]
        q[:3] = np.dot(A.T, t)[[2, 0, 1]]
        q[4] = 0
        q[3], q[5] = w
        return A

    def _to_ellipsoid_step2(self):
        q = self._q
        if (abs(q[3]) < self._prec and
            abs(q[5]) > self._prec):
            q[3] = q[5]
            q[5] = 0
            A = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]], dtype='double')
            t = q[[1, 2, 0]]
            q[:3] = np.dot(A.T, t)[[2, 0, 1]]
            return A
        else:
            return np.eye(3, dtype='double')
        
    def _to_ellipsoid_step3(self):
        q = self._q
        if (abs(q[3]) > self._prec and
            abs(q[1]) > self._prec):
            A = np.eye(3, dtype='double')
            A[0, 2] = -q[1] / (2 * q[3])
            q[0] -= (q[1] ** 2) / (4 * q[3])
            q[1] = 0
            return A
        else:
            return np.eye(3, dtype='double')

    def _to_ellipsoid_step3(self):
        q = self._q
        if (abs(q[5]) > self._prec and
            abs(q[2]) > self._prec):
            A = np.eye(3, dtype='double')
            A[1, 2] = -q[2] / (2 * q[5])
            q[0] -= (q[2] ** 2) / (4 * q[5])
            q[2] = 0
            return A
        else:
            return np.eye(3, dtype='double')
        
