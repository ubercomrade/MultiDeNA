from setuptools import setup
from Cython.Build import cythonize

setup(
    name='speedup',
    ext_modules=cythonize("lib.speedup.pyx", language_level = "3"),
    zip_safe=False,
)
