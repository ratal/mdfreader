from setuptools import setup, find_packages
import pkg_resources
from codecs import open  # To use a consistent encoding
from os import path
from distutils.extension import Extension
from distutils.version import LooseVersion
from warnings import warn

numpy_incl = pkg_resources.resource_filename('numpy', 'core/include')

# cython installed ?
min_cython_ver = '0.21'
try:
    import Cython
    ver = Cython.__version__
    use_cython = ver >= LooseVersion(min_cython_ver)
    from Cython.Build import cythonize
    ext = '.pyx' if use_cython else '.c'
    ext_modules = cythonize(Extension('dataRead', ['dataRead' + ext],
                                      include_dirs=[numpy_incl]))
except ImportError:
    use_cython = False

name = 'mdfreader'
version = '2.7.6'

description = 'A Measured Data Format file parser'

here = path.abspath(path.dirname(__file__))
# Get the long description from the relevant file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
long_description = long_description

# The project's main homepage.
url = 'https://github.com/ratal/mdfreader'

# Author details
author = 'Aymeric Rateau'
author_email = 'aymeric.rateau@gmail.com'

# Choose your license
license = 'GPL3'

# See https://pypi.python.org/pypi?%3Aaction=list_classifiers
classifiers = [
    # How mature is this project? Common values are
    #   3 - Alpha
    #   4 - Beta
    #   5 - Production/Stable
    'Development Status :: 4 - Beta',

    # Indicate who your project is intended for
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering',

    # Pick your license as you wish (should match "license" above)
    'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

    # Specify the Python versions you support here. In particular, ensure
    # that you indicate whether you support Python 2, Python 3 or both.
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6'
]

# What does your project relate to?
keywords = 'Parser MDF file'

# You can just specify the packages manually here if your project is
# simple. Or you can use find_packages().
packages = find_packages(exclude=['contrib', 'docs', 'tests*'])

# List run-time dependencies here.  These will be installed by pip when your
# project is installed. For an analysis of "install_requires" vs pip's
# requirements files see:
# https://packaging.python.org/en/latest/technical.html#install-requires-vs-requirements-files
install_requires = ['numpy>=1.14', 'sympy', 'lxml']

# List additional groups of dependencies here (e.g. development dependencies).
# You can install these using the following syntax, for example:
# $ pip install -e .[dev,test]
extras_require = {
    'export': ['scipy', 'h5py', 'xlwt', 'xlwt3', 'openpyxl>2.0', 'pandas'],
    'plot': ['matplotlib'],
    'converter': ['PyQt4'],
    'experimental': ['bitarray'],
    'compression': ['blosc'],
}

# If there are data files included in your packages that need to be
# installed, specify them here.  If using Python 2.6 or less, then these
# have to be included in MANIFEST.in as well.
# package_data={
#    'sample': ['package_data.dat'],
#},

# Although 'package_data' is the preferred approach, in some case you may
# need to place data files outside of your packages.
# see http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files
# In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
# data_files=[('my_data', ['data/data_file'])],

# To provide executable scripts, use entry points in preference to the
# "scripts" keyword. Entry points provide cross-platform support and allow
# pip to create the appropriate form of executable for the target platform.
entry_points = {
    'console_scripts': [
    'mdfconverter=mdfconverter.mdfconverter:main',
    ],
}

if use_cython:
    setup(name=name, version=version, description=description, long_description=long_description,
          url=url, author=author, author_email=author_email, license=license, classifiers=classifiers,
          keywords=keywords, packages=packages, install_requires=install_requires, extras_require=extras_require,
          entry_points=entry_points, ext_modules=ext_modules)
    warn('It is strongly advised to install Cython for performance and robustness purpose')
else:
    extras_require.pop('experimental')
    install_requires.append('bitarray')  # replaces cython requirement by bitarray
    setup(name=name, version=version, description=description, long_description=long_description,
          url=url, author=author, author_email=author_email, license=license, classifiers=classifiers,
          keywords=keywords, packages=packages, install_requires=install_requires, extras_require=extras_require,
          entry_points=entry_points)
