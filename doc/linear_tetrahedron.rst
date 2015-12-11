Linear tetrahedron method
==========================

This document is a summary of the paper,
MacDonald *et al.*, J. Phys. C: Solid State Phys., **12**, 2991
(1979).

For computer implementation of the tetrahedron method, the paper,
Peter E. Blöchl *et al*., Phys. Rev. B **49**, 16223 (1994), 
is really useful.

Functions
----------

A general spectral function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   G(\omega) = \int_\mathrm{BZ} \frac{d\mathbf{k}
   F(\mathbf{k})}{(\omega -\omega(\mathbf{k})-i\epsilon)}

Density of states :math:`g(\omega)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   g(\omega) = V_\mathrm{MZ} \sum^{N}_{i=1} g(\omega, \omega^i_1,
   \omega^i_2, \omega^i_3, \omega^i_4) \equiv V_\mathrm{MZ}
   \sum^{N}_{i=1} g^i

where :math:`N` is the number of tetrahedral microzones, :math:`i` is
the index running throughout the tetrahedral microzones, and
:math:`V_\mathrm{MZ}` is the microzone volume. 

Number of particles :math:`n(\omega)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   n(\omega) \equiv \int^\omega_{-\infty} d\omega'g(\omega')
   = V_\mathrm{MZ} \sum^{N}_{i=1} n(\omega, \omega^i_1,
   \omega^i_2, \omega^i_3, \omega^i_4) \equiv V_\mathrm{MZ}
   \sum^{N}_{i=1} n^i


.. _imagpart_func:
   
Imaginary part of a spectral function :math:`I(\omega)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   I(\omega) = \int d\mathbf{k} F(\mathbf{k})\delta(\omega
   -\omega(\mathbf{k})) = V_\mathrm{MZ} \sum^{N}_{i=1} g^i
   \sum^4_{k=1} I_k(\omega, \omega^i_1, \omega^i_2, \omega^i_3,
   \omega^i_4) F^i_k
   \equiv V_\mathrm{MZ} \sum^{N}_{i=1} g^i \sum^4_{k=1} I^i_k F^i_k

.. _integrated_func:

Integrated function of imaginary part of a spectral function :math:`J(\omega)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   J(\omega) = \int^\omega_{-\omega} d\omega' I(\omega') =
   \int_\mathrm{BZ} d\mathbf{k} F(\mathbf{k})\theta(\omega
   -\omega(\mathbf{k})) = V_\mathrm{MZ} \sum^{N}_{i=1} n^i
   \sum^4_{k=1} J_k(\omega, \omega^i_1, \omega^i_2, \omega^i_3,
   \omega^i_4) F^i_k \equiv V_\mathrm{MZ} \sum^{N}_{i=1} n^i
   \sum^4_{k=1} J^i_k F^i_k

Tetrahedron values
-------------------

:math:`\omega < \omega^i_1`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   g(\omega) &= 0 \\
   n(\omega) &= 0 \\
   I(\omega) &= 0 \\
   J(\omega) &= 0   

:math:`\omega^i_1 < \omega < \omega^i_2`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Vertices of a tetrahedron are :math:`\mathbf{k}_1`,
:math:`\mathbf{k}_\alpha = f_{21}\mathbf{k}_2 + f_{12}\mathbf{k}_1`,
:math:`\mathbf{k}_\beta = f_{31}\mathbf{k}_3 +  f_{13}\mathbf{k}_1`,
:math:`\mathbf{k}_\gamma = f_{41}\mathbf{k}_4 +  f_{14}\mathbf{k}_1`,
where :math:`f_{nm} \equiv (\omega - \omega_m) / (\omega_n - \omega_m)`.

.. math::

   \int_\mathrm{MZ} \tilde{F}(\mathbf{k}) \theta(\omega -
   \tilde{\omega}(\mathbf{k})) = V_\mathrm{MZ} n^i
   \tilde{F}((\mathbf{k}_1 + \mathbf{k}_\alpha + \mathbf{k}_\beta +
   \mathbf{k}_\gamma) / 4) = V_\mathrm{MZ} n^i \sum^4_{k=1} J^i_k
   F^i_k
   
.. math::

   n^i =& f_{21} f_{31} f_{41} \\
   J^i_1 =& (1 + f_{12} + f_{13} + f_{14})/4 \\
   J^i_k =& f_{k1} / 4\;\;\;\;\;\;\; k = 2, 3, 4 \\

By differentiating with respect to :math:`\omega`,

.. math::

   g^i =& 3n^i / (\omega - \omega_1) = 3 f_{21} f_{31} f_{41} /
   (\omega - \omega_1) = 3 f_{21} f_{31} / \Delta_{41} \\
   I^i_1 =& (f_{12} + f_{13} + f_{14})/3 \\
   I^i_k =& f_{k1} / 3\;\;\;\;\;\;\; k = 2, 3, 4 \\

where :math:`\Delta_{nm} \equiv \omega_n - \omega_m`.
   
:math:`\omega^i_2 < \omega < \omega^i_3`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The occupied part of the tetrahedron is the sum of the following three
tehtrahedra. Vertices of the first tetrahedron are
:math:`\mathbf{k}_1`, :math:`\mathbf{k}_2`, :math:`\mathbf{k}_\alpha =
f_{42}\mathbf{k}_4 + f_{24}\mathbf{k}_2`, :math:`\mathbf{k}_\beta =
f_{32}\mathbf{k}_3 + f_{23}\mathbf{k}_2`, and volume
:math:`V_{12\alpha\beta} = V_\mathrm{MZ}f_{42}f_{32}`. Vertices of the
second tetrahedron are :math:`\mathbf{k}_1`,
:math:`\mathbf{k}_\alpha`, :math:`\mathbf{k}_\beta`,
:math:`\mathbf{k}_\delta = f_{41}\mathbf{k}_4 + f_{14}\mathbf{k}_1`,
and volume :math:`V_{1\alpha\beta\delta} =
V_\mathrm{MZ}f_{41}f_{24}f_{32}`. Vertices of the third tetrahedron
are :math:`\mathbf{k}_1`, :math:`\mathbf{k}_\beta`,
:math:`\mathbf{k}_\delta`, :math:`\mathbf{k}_\gamma =
f_{31}\mathbf{k}_3 + f_{13}\mathbf{k}_1`, and volume
:math:`V_{1\beta\delta\gamma} = V_\mathrm{MZ}f_{41}f_{31}f_{23}`.

.. math::

   n^i =& f_{42} f_{32} + f_{41} f_{24} f_{32} + f_{41} f_{31} f_{23}
   = (V_{12\alpha\beta} + V_{1\alpha\beta\delta} +
   V_{1\beta\delta\gamma}) / V_\mathrm{MZ} \\
   J^i_1 =& (V_{12\alpha\beta} + V_{1\alpha\beta\delta}(1 + f_{14}) +
   V_{1\beta\delta\gamma} (1 + f_{14} + f_{13})) / 4 n^i V_\mathrm{MZ} \\
   J^i_2 =& (V_{12\alpha\beta} (1 + f_{24} + f_{23}) +
   V_{1\alpha\beta\delta}(f_{24} + f_{23}) +
   V_{1\beta\delta\gamma} f_{23}) / 4 n^i V_\mathrm{MZ} \\ 
   J^i_3 =& (V_{12\alpha\beta} f_{32} +
   V_{1\alpha\beta\delta} f_{32} +
   V_{1\beta\delta\gamma} (f_{32} + f_{31})) / 4 n^i V_\mathrm{MZ} \\ 
   J^i_4 =& (V_{12\alpha\beta} f_{42} +
   V_{1\alpha\beta\delta} (f_{42} + f_{41}) +
   V_{1\beta\delta\gamma} f_{41}) / 4 n^i V_\mathrm{MZ} \\
   g^i =& 3\Delta^{-1}_{41} (f_{23} f_{31} + f_{32} f_{24}) \\
   I^i_1 =& f_{14}/3 + f_{13} f_{31} f_{23} / g^i \Delta_{41} \\
   I^i_2 =& f_{23}/3 + f^2_{24} f_{32} / g^i \Delta_{41} \\
   I^i_3 =& f_{32}/3 + f^2_{31} f_{23} / g^i \Delta_{41} \\
   I^i_4 =& f_{41}/3 + f_{42} f_{24} f_{32} / g^i \Delta_{41}

:math:`\omega^i_3 < \omega < \omega^i_4`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The occupied part of the tetrahedron is the full tetrahedron minus the
tetrahedron with vertices at :math:`\mathbf{k}_4`,
:math:`\mathbf{k}_\beta = f_{42}\mathbf{k}_4 + f_{24}\mathbf{k}_2`,
:math:`\mathbf{k}_\delta = f_{41}\mathbf{k}_4 + f_{14}\mathbf{k}_1`
and :math:`\mathbf{k}_\gamma = f_{43}\mathbf{k}_4 +
f_{34}\mathbf{k}_3` and volume :math:`V_{4\beta\delta\gamma} =
V_\mathrm{MZ}f_{14}f_{24}f_{34}`.

.. math::

   n^i =& (1 - f_{14} f_{24} f_{34}) \\
   J^i_1 =& (1 - f^2_{14} f_{24} f_{34}) / 4n^i \\
   J^i_2 =& (1 - f_{14} f^2_{24} f_{34}) / 4n^i \\
   J^i_3 =& (1 + f_{14} f_{24} f^2_{34}) / 4n^i \\
   J^i_4 =& [1 - f_{14} f_{24} f_{34}(1 + f_{41} + f_{42} + f_{43})] / 4n^i \\
   g^i =& 3(1 - n^i)/(\omega_4 - \omega) = 3 f_{24} f_{34} / \Delta_{41} \\
   I^i_k =& f_{k4} / 3 \;\;\;\;\;\;\; k = 1, 2, 3 \\
   I^i_4 =& (f_{41} + f_{42} + f_{43}) / 3 \\

:math:`\omega^i_4 < \omega`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   n^i =& 1 \\
   J^i_k =& 1/4 \;\;\;\;\;\;\; k = 1, 2, 3, 4 \\

Implementation
---------------

..
   Each grid point on a uniform mesh (:math:`n_1 \times n_2 \times n_3`)
   is supposed to located at the origin of a parallelepiped microzone
   formed by reciprocal lattice vectors :math:`\mathbf{a}^*`,
   :math:`\mathbf{b}^*`, and :math:`\mathbf{c}^*`.

Following Blöchl's paper, we focus on a grid point :math:`p` and the
integral is made for the grid points. This is useful because the
symmetry handling is quite easy and it is treated in the same way as
the smearing method. Along this way, we rearrange the original
integral :math:`\sum^{N}_{i=1} g^i \sum^4_{k=1} I^i_kF^i_k` shown in
the section :ref:`imagpart_func` of this resume to
:math:`\sum^{N_\mathrm{mesh}}_{p=1}\sum^{24}_{l(p)=1} \sum^4_{k=1}
g^{l(p)} I^{l(p)}_k F^{l(p)}_k\delta_{kp}`, where for the original
integral, :math:`N` is the number of tetrahedral microzones, :math:`i`
is the index running throughout the tetrahedral microzones, and
:math:`k` gives the vertices of each tetrahedron. For the rearranged
integral, :math:`N_\mathrm{mesh}` is the number of grid points on the
uniform mesh, :math:`N = 6N_\mathrm{mesh}`, and :math:`l(p)` is the
tetrahedral microzone that runs through 24 tetrahedral microzones
chosen so that one of their vertices is located on our focused grid
point :math:`p` in each tetrahdron. In the computation, although the
index :math:`k` runs through four vertices of each tetrahedron, 
only one of them that correponds to the grid point :math:`p` is taken
to the sum.

1. Each parallelepiped microzone is divided into six tetrahedra by choosing the
   shortest main diagonal of the parallelepiped microzone as the common edge of
   all six tetrahedra.

2. A grid point :math:`p` is shared by 24 different tetrahedra, and
   one vertex for each tetrahedron. These 24 tetrahedra are chosen
   among 8 parallelepiped microzones depending on the choice of the
   main diagonal.

3. Among 48 tetrahedra, 24 tetrahedra are chosen by searching
   tetrahedra that contain the grid point :math:`p`. The 24 tetrahedra are
   stored in a (24, 4) array. In each row, index of the element of the
   grid point is remembered in a list of 24 elements.

4. Energies (or any real values) at grid points given by the (24, 4)
   tetrahedra array are input from outside and the energies of each
   row are stored in another (24, 4) array in ascending order. In each
   row of the energy array, correspondence between the grid point and
   its energy has to be kept as an index.

