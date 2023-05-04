from setuptools import Extension, setup

module = Extension("myspkmeans", sources=['spkmeansmodule.c' ,'wam.c', 'ddg.c', 'gl.c', 'jacobi.c', 'spkmeans.c'])
setup(name='myspkmeans',
     version='1.0',
     description='Python wrapper for custom C extension',
     ext_modules=[module])
