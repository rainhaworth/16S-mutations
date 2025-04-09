from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
    name='Faster gapped kmer search',
    ext_modules=cythonize("utils/find_gapped.pyx"),
    include_dirs=[numpy.get_include()]
)