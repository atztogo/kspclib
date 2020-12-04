import os
import numpy
from setuptools import setup, Extension

if os.path.exists('src'):
    source_dir = "src"
else:
    source_dir = os.path.join(os.pardir, "src")

include_dirs = [source_dir, numpy.get_include()]
sources = [os.path.join(source_dir, filename) for filename
           in ('kspclib.c',
               'kpoint.c',
               'tetrahedron_method.c',
               'kgengrid.c',
               'kgrid.c',
               'mathfunc.c')]
extra_compile_args = []
extra_link_args = []
define_macros = []

version_nums = [None, None, None]
with open(os.path.join(source_dir, "version.h")) as w:
    for line in w:
        for i, chars in enumerate(("MAJOR", "MINOR", "MICRO")):
            if chars in line:
                version_nums[i] = int(line.split()[2])

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
