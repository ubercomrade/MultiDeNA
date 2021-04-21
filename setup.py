"""pipeline for chip-seq analisys"""
from setuptools import setup, Extension
from Cython.Build import cythonize

with open("README.md", "r") as fh:
    long_description = fh.read()

install_requires = ['setuptools>=18.0', 'cython', 'numpy', 'scipy', 'matplotlib']
ext_modules = [Extension('lib.speedup', sources=['lib/speedup.pyx'],)]


#cmdclass={'build_ext': Cython.Build.build_ext},
setup(
    name='pipeline',
    version='0.0.1',
    description='pipeline for chip-seq analisys',
    author='Anton Tsukanov',
    author_email='tsukanov@bionet.nsc.ru',
    url='http://github.com/ubercomrade/pipeline',
    package_dir={'lib' : 'lib'},
    packages=[
        'lib',
        'tools'
    ],
    scripts=['pipeline.py',
    'scripts/cross-validation.py',],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        'Programming Language :: Cython',
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    zip_safe=False,
    install_requires=install_requires,
    setup_requires=install_requires,
    python_requires='>=3.6',
    ext_modules = cythonize(ext_modules)
)