from distutils.core import setup, Extension
from Cython.Distutils import build_ext
spglibroot="/home/flagre/localprogs/spglib-1.7.0"

ext = [Extension("spglib_interface", ["spglib_interfacec.pyx"],
                 include_dirs=[spglibroot+"/include/spglib"],
                 library_dirs=[spglibroot+"/lib/static"],
                 libraries=["symspg", "stdc++"],
		 language="c++"),
       Extension("tetrahedron_method_interface", ["kgrid.c","tetrahedron_method.c","tetrahedron_method_interfacec.pyx"],
                 libraries=["stdc++"],
		 language="c++")]

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=ext)
