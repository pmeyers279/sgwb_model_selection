#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) Patrick Meyers (2013)
#
# This file is part of SeisMon
#
# SeisMon is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SeisMon is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SeisMon.  If not, see <http://www.gnu.org/licenses/>

"""Setup script for SeisMon
"""

import glob
import os.path
from setuptools import (find_packages, setup)

from utils import version

PACKAGENAME = 'sbms'

VERSION_PY = os.path.join(PACKAGENAME, 'version.py')

# set version information
vcinfo = version.GitStatus()
vcinfo(VERSION_PY)

DESCRIPTION = 'SGWB Model Selection'
LONG_DESCRIPTION = ''
AUTHOR = 'Patrick Meyers'
AUTHOR_EMAIL = 'patrick.meyers@ligo.org'
LICENSE = 'GPLv3'

# VERSION should be PEP386 compatible (http://www.python.org/dev/peps/pep-0386)
VERSION = vcinfo.version

# Indicates if this version is a release version
RELEASE = vcinfo.version != vcinfo.id and 'dev' not in VERSION

# Use the find_packages tool to locate all packages and modules
packagenames = find_packages(exclude=['utils'])

# find all scripts
scripts = glob.glob('bin/*')
import unittest
def my_test_suite():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('sbms/tests', pattern='*test.py')
    return test_suite

setup(name=PACKAGENAME,
      version=VERSION,
      description=DESCRIPTION,
      scripts=scripts,
      packages=packagenames,
      ext_modules=[],
      requires=[],
      provides=[PACKAGENAME],
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      long_description=LONG_DESCRIPTION,
      zip_safe=False,
      test_suite='setup.my_test_suite',
      use_2to3=True
      )
