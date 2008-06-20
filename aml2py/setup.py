# import
from distutils.core import setup

LONG_DESC = """ aml2py convert AML files to python files."""

version= "0.1.2"

setup(name="aml2py",
      version=version,
      description="Convert aml files to python files",
      author="Christophe Pradal",
      author_email="christophe.pradal@cirad.fr",
      url="http://openalea.gforge.inria.fr",
      license="GPL",
      long_description=LONG_DESC,
      py_modules=['aml2py'])

