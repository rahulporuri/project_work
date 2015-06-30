from distutils.core import setup
from Cython.Build import cythonize

setup(
	name = 'routines',
	ext_modules = cythonize('routines.pxd'),
)
