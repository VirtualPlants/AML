num_jobs = 8



def complete_path(path):
	import os

	def find_version_dir(path):
		exts = os.listdir(path)
		def mycmp(a, b):
		   a = map(int, a.split('.'))
		   b = map(int, b.split('.'))
		   for ai, bi in zip(a, b):
		   	  if ai == bi: continue
		   	  return cmp(ai, bi)
		   return 0
		exts.sort(mycmp)
		return exts[0]
	prefix, suffix = path.split('/VERSION/')
	return os.path.join(prefix, find_version_dir(prefix), suffix)

# if qhull installed from brew, with default Qt :
qhull_libs_suffix = '.6'
qhull_lib = complete_path('/usr/local/Cellar/qhull/VERSION/lib')
qhull_includes = complete_path('/usr/local/Cellar/qhull/VERSION/include')

QTDIR = '/usr/local'
WITH_CGAL = True
WITH_MPFR = WITH_GMP = WITH_ANN = True
ann_includes = '/usr/local/include'
ann_libpath = '/usr/local/lib'
boost_includes = mpfr_includes = gmp_includes = '/usr/local/include'
boost_libpath = mpfr_libpath = gmp_libpath = '/usr/local/lib'
boost_libs_suffix = '-mt'

eigen_includes = complete_path('/usr/local/Cellar/eigen/VERSION/include/eigen3')

EXTRA_CCFLAGS = '-stdlib=libc++ -fno-inline -DWITH_QHULL_2011 -DCGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES -DCGAL_EIGEN3_ENABLED -I%s' % eigen_includes
EXTRA_LINK_FLAGS = '-F/usr/local/lib'

# If X11 has been installed manually :
gl_includes = '/usr/X11/include'

bison_bin = complete_path('/usr/local/Cellar/bison/VERSION/bin')
flex_bin = complete_path('/usr/local/Cellar/flex/VERSION/bin')

bison_libpath = complete_path('/usr/local/Cellar/bison/VERSION/lib')
flex_libpath = complete_path('/usr/local/Cellar/flex/VERSION/lib')

bison_include = complete_path('/usr/local/Cellar/bison/VERSION/include')
flex_include = complete_path('/usr/local/Cellar/flex/VERSION/include')

