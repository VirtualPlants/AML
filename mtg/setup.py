from setuptools import setup, find_packages
import os, sys
from os.path import join as pj
 
packagename = 'mtg'
namespace = 'openalea'
build_prefix = "build-scons"

# Scons build directory
scons_parameters=["build_prefix="+build_prefix]


# dependencies
#install_requires = ['vplants.amlobj', 'plantgl']
install_requires = []
setup_requires = ['openalea.deploy'] + install_requires

if __name__ == '__main__':
    
    setup(name='VPlants.Mtg',
          version='0.3.0',
          author='Christophe Godin',
          description='mtg structure',
          url='http://www-sop.inria.fr/virtualplants/',
          license='GPL',
 
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
          dependency_links = ['http://openalea.gforge.inria.fr/pi'],
          )
