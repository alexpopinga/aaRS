"""Build fast algorithms"""
from setuptools import setup, Extension
import numpy

setup(
    ext_modules=[
        Extension("aars_algorithms_fast", ["aars_algorithms_fast.pyx"],
                  include_dirs=[numpy.get_include()])
    ]
)
