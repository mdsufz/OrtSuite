#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Setup file for OrtSuite.

    This file was generated with PyScaffold 3.1.
    PyScaffold helps you to put up the scaffold of your new Python project.
    Learn more under: https://pyscaffold.org/
"""
import sys

from pkg_resources import require, VersionConflict
from setuptools import setup, find_packages

try:
    require('setuptools>=38.3')
except VersionConflict:
    print("Error: version of setuptools is too old (<38.3)!")
    sys.exit(1)

with open('README.md') as f:
    long_description = f.read()

if __name__ == "__main__":
    setup(use_pyscaffold=True,
          name='OrtSuite',
          version='0.01',
          description='Prediction of putative microbial interactions.',
          licence="GNU-LPL",  # Check license
          long_description=long_description,
          author='Microbial Systems Data Science group',
          author_email='msds@ufz.de',
          url='',
          # Commands
          scripts=['OrtAn/src/aux.py', 'OrtAn/src/diamond_mp.py','OrtScraper/src/MultipleRequests.py', 'OrtScraper/src/aux.py'],
          # classifiers=[],
          packages=find_packages(),
          install_requires=[],
          entry_points={
              'console_scripts':
                  ['create_project=OrtAn/src.create_project:run',
                   'delete_project=OrtAn/src.delete_project:run',
                   'relaxed_search=OrtAn/src.relaxed_search:run',
                   'restrictive_search=OrtAn/src.restrictive_search:run',
                   'check_projects=OrtAn/src.check_projects:run',
                   'annotation=OrtAn/src.annotation:run',
                   'create_db=OrtAn/src.create_db:run',
                   'download_kos=OrtScraper/src.download_kos:run'
                   ],
          },
          # Data - see this better
          # package_data=[]
          )
