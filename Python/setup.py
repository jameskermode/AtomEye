#!/usr/bin/env python
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2010
# HQ X
# HQ X   These portions of the source code are released under the GNU General
# HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HQ X
# HQ X   If you would like to license the source code under different terms,
# HQ X   please contact James Kermode, james.kermode@gmail.com
# HQ X
# HQ X   When using this software, please cite the following reference:
# HQ X
# HQ X   http://www.jrkermode.co.uk/quippy
# HQ X
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

import sys, os, re, string
from numpy.distutils.core import setup, Extension

major, minor = sys.version_info[0:2]
if (major, minor) < (2, 4):
    sys.stderr.write('Python 2.4 or later is needed to use this package\n')
    sys.exit(1)

atomeye_dir =  os.path.join(os.getcwd(), '..')
atomeye_libs = ['m', 'Xpm', 'Xext', 'X11', 'mpi', 'png', 'z', 'jpeg', 'history', 'ncurses',
                'netcdf', 'hdf5', 'hdf5_hl', 'readline', 'curl']
atomeye_libdirs = ['/opt/local/lib']
atomeye_extra_link_args = ['-framework Accelerate']

os.system('cd .. && make atomeyelib')

setup(name='atomeye',
      py_modules = ['atomeye'],
      ext_modules = [Extension(name='_atomeye',
                               sources=['atomeyemodule.c'],
                               library_dirs=atomeye_libdirs + [os.path.join(atomeye_dir, 'lib')],
                               libraries=['AtomEye', 'AX', 'Atoms', 'VecMat3', 'VecMat', 'IO', 'Scalar', 'Timer'] + atomeye_libs,
                               include_dirs=[os.path.join(atomeye_dir,'include')],
                               define_macros=[],
                               extra_link_args=atomeye_extra_link_args,
                               depends=[os.path.join(atomeye_dir, 'lib/libAtomEye.a')])
                     ],
      version=os.popen('svnversion -n .').read(),
      description='Python bindings to AtomEye',
      author='James Kermode',
      author_email='james.kermode@kcl.ac.uk',
      url='http://www.jrkermode.co.uk/AtomEye')
