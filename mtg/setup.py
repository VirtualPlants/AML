# -*- coding: utf-8 -*-
"""setup for mtg package"""
__revision__ = "$Id$"

import os, sys
from os.path import join as pj

from setuptools import setup
from openalea.deploy.binary_deps import binary_deps
#from openalea.deploy.metainfo import read_metainfo

# metadata = read_metainfo('metainfo.ini', verbose=True)
# for key,value in metadata.iteritems():
#     exec("%s = '%s'" % (key, value))

# Metainformation
version = "1.2.0"
name = "OpenAlea.amlmtg"
description="mtg structure"
authors="Christophe Godin"
url="https://github.com/VirtualPlants/AML/tree/master/mtg"

# Scons build directory
build_prefix = "build-scons"
scons_parameters = ["build_prefix="+build_prefix]


# dependencies
install_requires = [binary_deps('openalea.amlobj'), 
                    binary_deps('openalea.plantgl')]
install_requires = []
setup_requires = ['openalea.deploy'] 

if __name__ == '__main__':

    setup(name=name,
          version=version,
          author=authors,
          description=description,
          url=url,
          license=license,

          # Define where to execute scons
          scons_scripts=['SConstruct'],
          # Scons parameters  
          scons_parameters=scons_parameters,

          # Packages
          #packages=
          #package_dir={'':'build/lib'},

          # Add package platform libraries if any
          include_package_data=True,
          zip_safe = False,

          # Specific options of openalea.deploy
          lib_dirs = {'lib' : pj(build_prefix, 'lib'),},
          inc_dirs = { 'include' : pj(build_prefix, 'include') },
          

          # Dependencies
          setup_requires = setup_requires,
          install_requires = install_requires,
          #dependency_links = ['http://openalea.gforge.inria.fr/pi'],
          #pylint_packages  = ['.']
        )
