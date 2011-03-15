# -*- coding: utf-8 -*-
__revision__ = "$Id$"

# import
#from distutils.core import setup
from setuptools import setup, find_packages

from openalea.deploy.metainfo import read_metainfo

metadata = read_metainfo('metainfo.ini', verbose=True)
for key,value in metadata.iteritems():
    exec("%s = '%s'" % (key, value))



setup(name=name,
      version=version,
      description=description,
      long_description=long_description,
      author=author,
      author_email=author_email,
      url=url,
      license=license,

      py_modules=['aml2py'], 
      package_dir = {'':'src'}, 
      # entry_points
      entry_points = { 
        "console_scripts": ['aml2py = aml2py:main'],
        },
)

