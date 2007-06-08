# -+- python -+-

from openalea.distx import setup, find_packages, find_package_dir, Shortcut 
from os.path import join as pj
 
packagename = 'aml'
namespace = 'openalea'

# Scons build directory
build_prefix= "build-scons"
 
LONG_DESC = """
AmlPy provides a robust Python interface to the aml
programming language.  It can manage all kinds of AML objects and can
execute arbitrary AML functions. All the errors from the AML language are
converted to Python exceptions.
"""

if __name__ == '__main__':
    
    setup(name=packagename,
          version='2.0.0',
          author='',
          author_email='',
          description='aml and amlPylibraries',
          long_description=LONG_DESC,
          url='',
          license='GPL',
 
          # Define where to execute scons
          scons_scripts=['SConstruct'],
          # Scons parameters  
          scons_parameters=["build_prefix="+build_prefix],
 
        
          # Packages
          packages=find_packages(where='src', namespace=namespace),
          package_dir=find_package_dir(where='src', namespace=namespace), 
      
          # Add package platform libraries if any
          include_package_lib=True,
                    
	  # Copy shared data in default OpenAlea directory
	  # Map of 'destination subdirectory' : 'source subdirectory'
	  external_data={pj('include'):  pj(build_prefix, 'include'),
                         pj('lib'):  pj(build_prefix, 'lib'),
                         pj('bin'):  pj(build_prefix, 'bin'),
                         },
          )

