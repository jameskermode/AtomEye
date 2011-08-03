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
from distutils.sysconfig import parse_makefile

major, minor = sys.version_info[0:2]
if (major, minor) < (2, 4):
    sys.stderr.write('Python 2.4 or later is needed to use this package\n')
    sys.exit(1)

def expand_addsuffix(s):
    add_suffix = re.compile(r'\$[\(\{]addsuffix (.*?),(.*?)[\)\}]')
    try:
        m = add_suffix.search(s)
    except TypeError:
        return s
    if m:
        while m is not None:
            suffix, files = m.groups()
            s = add_suffix.sub(' '.join([f + suffix for f in files.split()]),s,1).strip()
            m = add_suffix.search(s)
    return s

def expand_addprefix(s):
    add_prefix =  re.compile(r'\$[\(\{]addprefix (.*?),(.*?)[\}\)]')
    try:
        m = add_prefix.search(s)
    except TypeError:
        return s
    if m:
        while m is not None:
            prefix, files = m.groups()
            s = add_prefix.sub(' '.join([prefix + f for f in files.split()]),s,1).strip()
            m = add_prefix.search(s)
    return s

atomeye_dir =  os.path.join(os.getcwd(), '..')
os.system('cd %s && make -p atomeyelib > %s/make.log' % (atomeye_dir, os.getcwd()))
makefile = parse_makefile('make.log')
os.unlink('make.log')

syslibs = makefile['SYSLIBS']
syslibs = expand_addprefix(syslibs)
syslibs = expand_addsuffix(syslibs)
syslibs = syslibs.split()

atomeye_libs = ['m', 'Xpm', 'Xext', 'X11', 'mpi', 'png', 'z', 'jpeg', 'history', 'ncurses', 'readline'] + [ f[2:] for f in syslibs if f.startswith('-l') ]
atomeye_libdirs = [ f[2:] for f in syslibs if f.startswith('-L') ]
atomeye_internal_libs = ['%s/lib%s.a' % (os.path.join(atomeye_dir, 'lib'), lib) for lib in [ 'AtomEye', 'AX', 'Atoms', 'VecMat3', 'VecMat', 'IO', 'Scalar', 'Timer'] ]
atomeye_extra_link_args = atomeye_internal_libs + [ f for f in syslibs if not f.startswith('-l') and not f.startswith('-L')]

quip_root_dir = os.environ['QUIP_ROOT']
quip_arch = os.environ['QUIP_ARCH']

if 'QUIPPY_LDFLAGS' in makefile:
    atomeye_extra_link_args.extend(makefile['QUIPPY_LDFLAGS'].split())

setup(name='atomeye',
      py_modules = ['atomeye'],
      ext_modules = [Extension(name='_atomeye',
                               sources=['atomeyemodule.c'],
                               library_dirs=atomeye_libdirs,
                               libraries=atomeye_libs,
                               include_dirs=[os.path.join(atomeye_dir,'include'), os.path.join(quip_root_dir,'libAtoms')],
                               define_macros=[],
                               extra_link_args=atomeye_extra_link_args,
                               depends=[os.path.join(atomeye_dir, 'lib/libAtomEye.a')])
                     ],
      version=os.popen('svnversion -n .').read(),
      description='Python bindings to AtomEye',
      author='James Kermode',
      author_email='james.kermode@kcl.ac.uk',
      url='http://www.jrkermode.co.uk/quippy/atomeye.html')
