from setuptools import setup, find_packages
from openalea.deploy.binary_deps import binary_deps
import os, sys
from os.path import join as pj
 
packagename = 'aml'
namespace = 'openalea'
build_prefix = "build-scons"
version = '0.7.0'

# Scons build directory
scons_parameters=["build_prefix="+build_prefix]


# platform dependencies
install_requires = [binary_deps('vplants.amlobj'), 
                    binary_deps('vplants.sequence_analysis'),
                    binary_deps('vplants.tree_matching'), 
                    binary_deps('vplants.mtg'), 
                    ]

if sys.platform.startswith('win'):
    install_requires += [binary_deps("boostpython"),"gnuplot"]

setup_requires = install_requires + ['openalea.deploy']

if __name__ == '__main__':
    
    setup(name='VPlants.Aml',
          version=version,
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

          packages=[namespace+".aml", namespace+".aml_wralea", 
                    namespace+".aml.amldoc",
                    namespace+".aml_wralea.wralea", 
                    namespace+".aml_wralea.demo", 
                    namespace+".aml_wralea.demo.change_point", 
                    namespace+".aml_wralea.demo.mtg", 
                    "vplants", "vplants.aml",
                    "amlPy"],

          package_dir={namespace+".aml" : "src/aml",
                       namespace+".aml_wralea" : "src/aml_wralea",
                       namespace+".aml.amldoc" : "src/aml/amldoc",
                       namespace+".aml_wralea.wralea" : "src/aml_wralea",
                       namespace+".aml_wralea.demo" : "src/aml_wralea/demo",
                       namespace+".aml_wralea.demo.change_point" : "src/aml_wralea/demo/change_point", 
                       namespace+".aml_wralea.demo.mtg": "src/aml_wralea/demo/mtg", 

                       "amlPy" : "src/amlPy",
                       "vplants": "src/vplants",
                       "" : 'src',
                       },
      
          # Add package platform libraries if any
          package_data = {namespace+".aml" : ['*.pyd', '*.so', '*.dylib'],},
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
            "wralea": ["aml = openalea.aml_wralea.wralea",
                       "demo = openalea.aml_wralea.demo"]
            },

          pylint_packages = ['.']
          )


    
