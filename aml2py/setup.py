# import
#from distutils.core import setup
from setuptools import setup, find_packages

LONG_DESC = """ aml2py convert AML files to python files."""

version= "0.7.0"

setup(name="aml2py",
      version=version,
      description="Convert aml files to python files",
      author="Christophe Pradal",
      author_email="christophe.pradal@cirad.fr",
      url="http://openalea.gforge.inria.fr",
      license="CECILL-C",
      long_description=LONG_DESC,
      py_modules=['aml2py'], 

      # entry_points
      entry_points = { 
        "console_scripts": ['aml2py = aml2py:main'],
        },
)

