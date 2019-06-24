from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

ext_modules = [Extension("make_pentaCY", ["penta_fncs.pyx"])]

setup(
	name = 'penta channel app',
	cmdclass = {'build_ext': build_ext},
	ext_modules = ext_modules
)
