from setuptools import setup, find_packages
import os, sys
from os.path import join as pj
 
packagename = 'aml'
namespace = 'openalea'
build_prefix = "build-scons"

# Scons build directory
scons_parameters=["build_prefix="+build_prefix]


# platform dependencies
if("win" in sys.platform):
    install_requires = ["boostpython","gnuplot"]
    setup_requires = install_requires 
else:
    install_requires = []
    setup_requires = []

# common deps
common_deps = ['vplants.amlobj', 'vplants.sequence_analysis',
         'vplants.tree_matching', 'vplants.flux', 'vplants.msvoxel',
         ]

if __name__ == '__main__':
    
    setup(name='VPlants.Aml',
          version='2.0.1',
          author='C. Godin, Y. Guedon, C.Pradal',
          description='aml and amlPy libraries',
          url='',
          license='GPL',
          
          # Define where to execute scons
          scons_scripts=['SConstruct'],
          # Scons parameters  
          scons_parameters=scons_parameters,
          
        
          # Packages
          namespace_packages = [namespace],
          create_namespaces = True,

          packages=[namespace+".aml", namespace+".aml.amldoc",namespace+".aml.wralea", "amlPy"],
          package_dir={namespace+".aml" : "src/aml",
                       namespace+".aml.amldoc" : "src/aml/amldoc",
                       namespace+".aml.wralea" : "src/aml/wralea",
                       "amlPy" : "src/amlPy"},
      
          # Add package platform libraries if any
          package_data = {namespace+".aml" : ['*.pyd', '*.so'],},
          include_package_data=True,
          zip_safe = False,

          # Specific options of openalea.deploy
          lib_dirs = {'lib' : pj(build_prefix, 'lib'),},
          inc_dirs = { 'include' : pj(build_prefix, 'include') },
          bin_dirs = { 'bin':  pj(build_prefix, 'bin') },
          

          # Dependencies
          setup_requires = setup_requires + common_deps + ['openalea.deploy'],
          install_requires = install_requires + common_deps,
          dependency_links = ['http://openalea.gforge.inria.fr/pi'],
          )


    
