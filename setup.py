from setuptools import setup
#from setuptools.extension import Extension

import numpy as np

import os
import sys
import subprocess

#from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
#from Cython.Distutils import build_ext


args = sys.argv[1:]

# Make a `cleanall` rule to get rid of intermediate and library files
if "cleanall" in args:
    print("Deleting cython files...")
    # Just in case the build directory was created by accident,
    # note that shell=True should be OK here because the command is constant.
    subprocess.Popen("rm -rf onezone_cython.egg-info",shell=True,executable="/bin/bash")
    subprocess.Popen("rm -rf onezone.egg-info",shell=True,executable="/bin/bash")
    subprocess.Popen("rm -rf dist",shell=True,executable="/bin/bash")

    for p in ['./','./onezone/cython_ext/']:
        subprocess.Popen("rm -rf " + p + "build", shell=True, executable="/bin/bash")
        subprocess.Popen("rm -rf " + p + "*.c", shell=True, executable="/bin/bash")
        subprocess.Popen("rm -rf " + p + "*.so", shell=True, executable="/bin/bash")


    sys.exit()

cython_extensions = [
    Extension(
          "onezone.cython_ext.sample_imf",
          ["onezone/cython_ext/sample_imf.pyx"],
        ),
    Extension(
           "onezone.cython_ext.cython_star",
          ["onezone/cython_ext/cython_star.pyx"],
        )
]

setup(
      name="onezone",
      version="0.1",
      author="Andrew Emerick",
      author_email="aemerick11@gmail.com",
      url="https://github.com/aemerick/onezone",
      packages = [
                  'onezone',
                  'onezone.analysis',
                  'onezone.plots',
                  'onezone.cython_ext'],
      install_requires=[
                'numpy',
                'scipy',
                'cython',
                'matplotlib',
                'h5py'], # and others ... need to finish
      ext_modules = cythonize(cython_extensions) # cythonize("onezone/cython_ext/sample_imf.pyx"),
#cython_extensions
)
