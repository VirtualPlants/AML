from setuptools import setup, find_packages
from openalea.deploy.binary_deps import binary_deps
import os, sys
from os.path import join as pj
 
packagename = 'aml'
namespace = 'openalea'
build_prefix = "build-scons"

# Scons build directory
scons_parameters=["build_prefix="+build_prefix]


# platform dependencies
install_requires = [binary_deps('vplants.amlobj'), 
                    binary_deps('vplants.sequence_analysis'),
                    binary_deps('vplants.tree_matching'), 
                    ]

if("win" in sys.platform):
    install_requires += [binary_deps("boostpython"),"gnuplot"]

setup_requires = install_requires + ['openalea.deploy']

if __name__ == '__main__':
    
    setup(name='VPlants.Aml',
          version='0.4.0',
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

          packages=[namespace+".aml", 
                    namespace+".aml.amldoc",
                    namespace+".aml.wralea", 
                    namespace+".aml.demo", 
                    namespace+".aml.demo.change_point", 
                    namespace+".aml.demo.mtg", 
                    "vplants.aml",
                    "amlPy"],

          package_dir={namespace+".aml" : "src/aml",
                       namespace+".aml.amldoc" : "src/aml/amldoc",
                       namespace+".aml.wralea" : "src/aml/wralea",
                       namespace+".aml.demo" : "src/aml/demo",
                       namespace+".aml.demo.change_point" : "src/aml/demo/change_point", 
                       namespace+".aml.demo.mtg": "src/aml/demo/mtg", 

                       "amlPy" : "src/amlPy",
                       "vplants.aml": "src/vplants/aml",
                       "" : 'src',
                       },
      
          # Add package platform libraries if any
          package_data = {namespace+".aml" : ['*.pyd', '*.so'],},
          include_package_data=True,
          zip_safe = False,

          # Specific options of openalea.deploy
          lib_dirs = {'lib' : pj(build_prefix, 'lib'),},
          inc_dirs = { 'include' : pj(build_prefix, 'include') },
          bin_dirs = { 'bin':  pj(build_prefix, 'bin') },
          share_dirs = { 'databases' : 'databases' },

          # Dependencies
          setup_requires = setup_requires ,
          install_requires = install_requires,
          dependency_links = ['http://openalea.gforge.inria.fr/pi'],

          # entry_points
          entry_points = {
            "wralea": ["aml = openalea.aml.wralea",
                       "demo = openalea.aml.demo"]
            },


          )


    
