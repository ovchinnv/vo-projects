from distutils.core import setup
from distutils.extension import Extension
import os
import sys
import platform

openmm_dir = '@OPENMM_DIR@'
openmmstruna_header_dir = '@OPENMMSTRUNA_HEADER_DIR@'
openmmstruna_library_dir = '@OPENMMSTRUNA_LIBRARY_DIR@'

# setup extra compile and link arguments on Mac
extra_compile_args = []
extra_link_args = []

if platform.system() == 'Darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl', '-rpath', openmm_dir+'/lib']

extension = Extension(name='_openmmstruna',
                      sources=['StrunaPluginWrapper.cpp'],
                      libraries=['OpenMM', 'OpenMMStruna'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), openmmstruna_header_dir],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), openmmstruna_library_dir],
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args
                     )

setup(name='OpenMMStruna',
      version='1.0',
      py_modules=['openmmstruna'],
      ext_modules=[extension],
     )
