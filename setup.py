"""pipeline for chip-seq analisys"""
from setuptools import setup, Extension
import Cython.Build

with open("README.md", "r") as fh:
    long_description = fh.read()

# Dependencies for pip:
install_requires = [
        'matplotlib >= 1.3.1',
        'numpy >= 1.18.5',
        'pandas >= 1.0.4',
        'cython >= 0.29.19'
]

# cython modules
ext_modules = [Extension(
        'lib.speedup',
        sources=['lib/speedup.pyx'],
    )]

setup(
    name='pipeline',
    version='0.0.1',
    description='pipeline for chip-seq analisys',
    author='Anton Tsukanov',
    author_email='tsukanov@bionet.nsc.ru',
    url='http://github.com/ubercomrade/pipeline',
    package_dir={'lib' : 'lib'},
    packages=[
        'lib.speedup',
        'tools'
    ],
    scripts=['pipeline.py'],
    cmdclass={'build_ext': Cython.Build.build_ext},
    install_requires=install_requires,
    setup_requires=install_requires,
    ext_modules=ext_modules,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        'Programming Language :: Cython',
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    zip_safe=False
)