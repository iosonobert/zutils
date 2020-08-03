#!/usr/bin/env python

from setuptools import setup, find_packages
from setuptools.dist import Distribution

class BinaryDistribution(Distribution):
    def is_pure(self):
        return False

		
setup(name='crude_readers',
      version='1.0.0',
      description='Package testing for 02',
      author='Andrew Zulberti',
      author_email='andrew.zulberti@gmail.com',
      #packages=['crude_readers'],
      packages=find_packages(),
      install_requires=['numpy','scipy','matplotlib','gsw','netcdf4','xarray'],
      license='unlicensed to all but author',
      include_package_data=True,
      distclass=BinaryDistribution,
    )
