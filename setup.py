"""Build script: all project metadata lives in pyproject.toml.

Only the optional Cython extension ``mdfreader.dataRead`` is declared here.
If it fails to compile, the package is still installed and falls back to the
slower pure-python code path (which needs ``bitarray``, see the ``fallback`` extra).
"""
from setuptools import setup, Extension

try:
    import numpy
    from Cython.Build import cythonize
    ext_modules = cythonize(
        [Extension("mdfreader.dataRead", ["mdfreader/dataRead.pyx"],
                   include_dirs=[numpy.get_include()],
                   optional=True)],  # a compilation failure does not abort the install
        language_level=3)
except ImportError:
    print("WARNING: Cython or numpy not available, dataRead extension not built. "
          "Install Cython with a compiler for performance and robustness.")
    ext_modules = []

setup(ext_modules=ext_modules)
