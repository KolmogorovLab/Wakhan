import os
import sys
import subprocess
import shutil

try:
    import setuptools
except ImportError:
    sys.exit("setuptools package not found. "
             "Please use 'pip install setuptools' first")

from setuptools import setup

# Make sure we're running from the setup.py directory.
script_dir = os.path.dirname(os.path.realpath(__file__))
if script_dir != os.getcwd():
    os.chdir(script_dir)

from src.__version__ import __version__

setup(name='wakhan',
      version=__version__,
      description='A tool for haplotype-specific somatic copy number aberrations/profiling from long reads sequencing data',
      url='https://github.com/KolmogorovLab/Wakhan',
      author='Tanveer Ahmad',
      author_email = 'tanveer.ahmad@nih.gov',
      license='MIT',
      packages=['src'],
      package_data={'src': ['annotations/*']},
      entry_points={'console_scripts': ['wakhan = src.main:main']},
      )