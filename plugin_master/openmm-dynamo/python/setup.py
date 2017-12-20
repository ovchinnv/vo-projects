from distutils.core import setup
from distutils.extension import Extension
import os
import sys
import platform

openmm_dir = '@OPENMM_DIR@'
openmmdynamo_header_dir = '@OPENMMDYNAMO_HEADER_DIR@'
openmmdynamo_library_dir = '@OPENMMDYNAMO_LIBRARY_DIR@'

# setup extra compile and link arguments on Mac
extra_compile_args = []
extra_link_args = []

if platform.system() == 'Darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl', '-rpath', openmm_dir+'/lib']

extension = Extension(name='_openmmdynamo',
                      sources=['DynamoPluginWrapper.cpp'],
                      libraries=['OpenMM', 'OpenMMDynamo'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), openmmdynamo_header_dir],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), openmmdynamo_library_dir],
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args
                     )

setup(name='OpenMMDynamo',
      version='1.0',
      py_modules=['openmmdynamo'],
      ext_modules=[extension],
     )
