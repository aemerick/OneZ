from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize

cython_extensions = [
    Extension(
        "cython_ext.sample_imf",
        ["cython_ext/sample_imf.pyx"],
        ),
]

setup(
      name="onezone",
      version="",
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
      ext_modules = cythonize("onezone/cython_ext/sample_imf.pyx"),
#cython_extensions
)
