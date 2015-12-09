from distutils.core import setup, Extension
from Cython.Distutils import build_ext

ext = [Extension("tetrahedron_method_interface", ["tetrahedron_method_interfacec.pyx"],
                           libraries=["tetrahedron", "stdc++"],
			   language="c++"),
setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=ext)
