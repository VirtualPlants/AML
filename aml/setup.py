# -+- python -+-

from openalea.distx import setup, find_packages, find_package_dir, Shortcut 
from os.path import join as pj
 
packagename = 'aml'
namespace = 'openalea'

# Scons build directory
build_prefix= "build-scons"
 

if __name__ == '__main__':
    
    setup(name=packagename,
          version='2.0.0',
          author='C. Godin, Y. Guédon, C.Pradal',
          author_email='',
          description='aml and amlPy libraries',
          url='',
          license='GPL',
 
          # Define where to execute scons
          scons_scripts=['SConstruct'],
          # Scons parameters  
          scons_parameters=["build_prefix="+build_prefix],
 
        
          # Packages
          packages=[namespace+".aml", "amlPy"],
          package_dir={namespace+".aml" : "src/aml",
                       "amlPy" : "src/amlPy"},
      
          # Add package platform libraries if any
          include_package_lib=True,
                    
	  # Copy shared data in default OpenAlea directory
	  # Map of 'destination subdirectory' : 'source subdirectory'
	  external_data={pj('include'):  pj(build_prefix, 'include'),
                         pj('lib'):  pj(build_prefix, 'lib'),
                         pj('bin'):  pj(build_prefix, 'bin'),
                         },
          )

