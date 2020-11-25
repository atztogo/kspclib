import numpy
from setuptools import setup, Extension

include_dirs = ['../src', numpy.get_include()]
sources = ['../src/kspclib.c',
           '../src/kpoint.c',
           '../src/tetrahedron_method.c',
           '../src/kgengrid.c',
           '../src/kgrid.c',
           '../src/mathfunc.c']
extra_compile_args = []
extra_link_args = []
define_macros = []

extension = Extension('kspclib._kspclib',
                      include_dirs=include_dirs,
                      sources=['_kspclib.c'] + sources,
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args,
                      define_macros=define_macros)

setup(name='kspclib',
      version='0.1',
      setup_requires=['numpy', 'setuptools'],
      description='This is the Kspclib module.',
      author='Atsushi Togo',
      author_email='atz.togo@gmail.com',
      url='https://github.com/atztogo/kspclib',
      packages=['kspclib'],
      install_requires=['numpy', ],
      provides=['kspclib'],
      platforms=['all'],
      ext_modules=[extension])
