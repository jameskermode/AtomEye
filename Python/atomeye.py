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

"""This module provides a high-level interface to the AtomEye extension module :mod:`_atomeye`. """

import sys
import numpy
import time
import imp
import os
import os.path
import shutil
import sys
import tempfile
from math import ceil, log10


def _tmp_pkg(dir=None):
    """
    Create a temporary package.

    Returns (name, path)
    """
    while True:
        path = tempfile.mkdtemp(dir=dir)
        name = os.path.basename(path)
        try:
            modinfo = imp.find_module(name)
            # if name is found, delete and try again
            os.rmdir(path)
        except:
            break
    init = file(os.path.join(path, '__init__.py'), 'w')
    init.close()
    return name, path


class MExt(object):
    """
    Load a unique copy of a module that can be treated as a "class instance".

    Adapted from code posted by Tod. A. Smith <Tod.A.Smith@aadc.com> to
    http://cens.ioc.ee/pipermail/f2py-users/2004-September/000921.html
    """

    def __init__(self, name):
        self.name = name
        # first find the "real" module on the "real" syspath
        srcfile, srcpath, srcdesc = imp.find_module(name)
        # now create a temp directory for the bogus package
        self._pkgname, self._pkgdir = _tmp_pkg('.')
        # copy the original module to the new package
        shutil.copy(srcpath, self._pkgdir)
        # import the module
        # __import__ returns the package, not the sub-module
        self._pkg = __import__(self._pkgname, globals(), locals(), [self.name])
        # return the module object
        self._module = getattr(self._pkg, self.name)
        # now add the module's stuff to this class
        self.__dict__.update(self._module.__dict__)

    def __del__(self):
        # remove module
        import sys, shutil
        del sys.modules[self._module.__name__]
        del sys.modules[self._pkg.__name__]
        # now try to delete the files and directory
        shutil.rmtree(self._pkgdir)
        # make sure the original module is loaded -
        # otherwise python crashes on exit
        # if MExt objects have not been explicitly 'del'd,
        # and __del__ is occurring on python shutdown, the import will fail
        # and the exception is caught here
        try:
            __import__(self.name)
        except:
            pass


default_settings = {'n->xtal_mode': 1,
                    'n->suppress_printout': 1,
                    'n->bond_mode': 1,
                    'n->atom_r_ratio': 0.5,
                    'key->BackSpace': 'load_config_backward'
                    }

name_map = {'positions': 'pos',
            'masses'   : 'mass',
            'numbers'  : 'Z' }

def on_click(iw, idx):
    if not iw in views:
        raise RuntimeError('Unexpected window id %d' % iw)
    views[iw].on_click(idx)
    

def on_advance(iw, mode):
    if not iw in views:
        raise RuntimeError('Unexpected window id %d' % iw)
    views[iw].on_advance(mode)

def on_close(iw):
    if not iw in views:
        raise RuntimeError('Unexpected window id %d' % iw)
    views[iw].on_close()


def on_new_window(iw):
    if iw in views:
        views[iw].is_alive = True
    else:
        views[iw] = AtomEyeView(window_id=iw)
    
class AtomEyeView(object):
    def __init__(self, atoms=None, window_id=None, copy=None, frame=0, delta=1, property=None, arrows=None, nowindow=False,
                 echo=False, block=False, verbose=True, *arrowargs, **arrowkwargs):
        self.atoms = atoms
        self.frame = frame
        self.delta = delta
        self.echo = echo
        self.block = block
        self.fortran_indexing = False
        self.verbose = verbose
        self.is_alive = False
        self.selection = []

        if window_id is None:
            self.start(copy, nowindow)
        else:
            self._window_id = window_id
            self.is_alive = True
            views[self._window_id] = self

        global view
        if view is None:
            view = self

        if property is not None or arrows is not None:
            self.show(obj=None, property=property, arrows=arrows, *arrowargs, **arrowkwargs)

    def _convert_atoms(self):
        self.current_atoms = None
        if self.atoms is None:
            title = '(null)'
            n_atom = 0
            cell = None
            arrays = None
        else:
            if hasattr(self.atoms, '__iter__'):
                self.current_atoms = self.atoms[self.frame]
                fmt = "%%0%dd" % ceil(log10(len(self.atoms)+1))
                title = 'AtomsList[%s] len=%s' % (fmt % self.frame, fmt % len(self.atoms))
            else:
                title = 'Atoms'
                self.current_atoms = self.atoms

            n_atom = len(self.current_atoms)

            self.fortran_indexing = False
            try:
                self.fortran_indexing = self.current_atoms.fortran_indexing
            except AttributeError:
                pass

            cell = self.current_atoms.get_cell()
            pbc = self.current_atoms.get_pbc()
            pos = self.current_atoms.positions
            
            for i, p in enumerate(pbc):
                if not p:
                    cell[i,i] = max(1.0, 2*(pos[:,i].max()-pos[:,i].min()))

            print 'Number of atoms: %d' % n_atom
            print 'Fortran indexing: %r' % self.fortran_indexing
            print 'Effective unit cell:'
            print cell

            try:
                arrays = self.current_atoms.properties
            except AttributeError:
                arrays = {}
                for key,value in self.current_atoms.arrays.iteritems():
                    arrays[name_map.get(key,key)] = value
            
        return title, n_atom, cell, arrays
        
            
    def start(self, copy=None, nowindow=False):
        if self.is_alive: return
        title, n_atom, cell, arrays = self._convert_atoms()

        icopy = -1
        if copy is not None:
            if isinstance(copy, AtomEye):
                icopy = copy._window_id
            elif isinstance(copy, int):
                icopy = copy
            else:
                raise TypeError('copy should be either an int or an AtomEye instance')

        # create our own unique version of the _atomeye extension module
        self._atomeye = MExt('_atomeye')
        self._atomeye.set_handlers(on_click, on_close, on_advance, on_new_window)

        self.is_alive = False
        self._window_id = self._atomeye.open_window(icopy,n_atom,cell,arrays,nowindow)

        views[self._window_id] = self
        while not self.is_alive:
            time.sleep(0.1)
        time.sleep(0.3)
        self._atomeye.set_title(self._window_id, title)
        self.update(default_settings)

    def on_click(self, idx):
        if self.current_atoms is None: return
        if idx >= len(self.current_atoms):
            idx = idx % len(self.current_atoms)
        if self.fortran_indexing:
            idx = idx + 1 # atomeye uses zero based indices

        self.selection.append(idx)

        if self.verbose:
            print
            try:
                self.current_atoms.print_atom(idx)
            except AttributeError:
                print self.current_atoms[idx]
            sys.stdout.flush()

    def on_advance(self, mode):
        if not hasattr(self.atoms,'__iter__'): return

        if mode == 'forward':
            self.frame += self.delta
        elif mode == 'backward':
            self.frame -= self.delta
        elif mode == 'first':
            self.frame = 0
        elif mode == 'last':
            self.frame = len(self.atoms)-1

        if self.frame > len(self.atoms)-1:
            try:
                self.atoms[self.frame]
            except IndexError:
                self.frame = self.frame % len(self.atoms)
                
        if self.frame < 1:
            self.frame = self.frame % len(self.atoms)

        print 'setting frame to %d' % self.frame
        sys.stdout.flush()
        self.show()

        if self.verbose:
            print
            print self.atoms[self.frame].params
            sys.stdout.flush()


    def on_close(self):
        global view
        
        self.is_alive = False
        if view is self:
            view = None
        del views[self._window_id]

    def show(self, obj=None, property=None, frame=None, arrows=None, *arrowargs, **arrowkwargs):
        if not self.is_alive:
            raise RuntimeError('is_alive is False')

        if obj is not None:
            self.atoms = obj
        
        if hasattr(self.atoms,'__iter__'):
            if frame is not None:
                if frame < 0: frame = len(self.atoms)-frame
                if frame >= len(self.atoms):
                    try:
                        self.atoms[self.frame]
                    except IndexError:
                        frame=len(self.atoms)-1
                self.frame = frame
       
        if property is not None:
            if isinstance(property,str):
                pass
            elif isinstance(property,int):
                if self.atoms.has_property('_show'):
                    self.atoms.remove_property('_show')
                self.atoms.add_property('_show', False)
                self.atoms._show[:] = [i == property for i in frange(self.atoms.n)]
                property = '_show'
            else:
                if self.atoms.has_property('_show'):
                    self.atoms.remove_property('_show')
                self.atoms.add_property('_show', property)
                property = '_show'

        title, n_atom, cell, arrays = self._convert_atoms()
        self._atomeye.load_atoms(self._window_id, title, n_atom, cell, arrays)
        if property is not None:
            self.aux_property_coloring(property)
        if arrows is not None:
            self.draw_arrows(arrows, *arrowargs, **arrowkwargs)

    def run_command(self, command):
        if not self.is_alive: 
            raise RuntimeError('is_alive is False')
        if self.echo:
            print command.strip()
        self._atomeye.run_command(self._window_id, command)
        if self.block:
            self.wait()

    def run_script(self, script):
        if type(script) == type(''):
            script = open(script)
            
        for line in script:
            self.run_command(line)
            self.wait()

    def __call__(self, command):
        self.run_command(command)

    def close(self):
        self.run_command('close')

    def update(self, D):
        for k, v in D.iteritems():
            self.run_command("set %s %s" % (str(k), str(v)))

    def save(self, filename):
        self.run_command("save %s" % str(filename))

    def load_script(self, filename):
        self.run_command("load_script %s" % str(filename))

    def key(self, key):
        self.run_command("key %s" % key)

    def toggle_coordination_coloring(self):
        self.run_command("toggle_coordination_coloring")

    def translate(self, axis, delta):
        self.run_command("translate %d %f " % (axis, delta))

    def shift_xtal(self, axis, delta):
        self.run_command("shift_xtal %d %f" % (axis, delta))

    def rotate(self, axis, theta):
        self.run_command("rotate %d %f" % (axis, theta))

    def advance(self, delta):
        self.run_command("advance %f" % delta)

    def shift_cutting_plane(self, delta):
        self.run_command("shift_cutting_plane %f" % delta)

    def change_bgcolor(self, color):
        self.run_command("change_bgcolor %f %f %f" % (color[0], color[1], color[2]))

    def change_atom_r_ratio(self, delta):
        self.run_command("change_atom_r_ratio %f" % delta)

    def change_bond_radius(self, delta):
        self.run_command("change_bond_radius %f" % delta)

    def change_view_angle_amplification(self, delta):
        self.run_command("change_view_angle_amplification %f" % delta)

    def toggle_parallel_projection(self):
        self.run_command("toggle_parallel_projection")

    def toggle_bond_mode(self):
        self.run_command("toggle_bond_mode" )

    def toggle_small_cell_mode(self):
        self.run_command("toggle_small_cell_mode")
        self.show()

    def normal_coloring(self):
        self.run_command("normal_coloring")

    def aux_property_coloring(self, auxprop):
        self.run_command("aux_property_coloring %s" % str(auxprop))

    def central_symmetry_coloring(self):
        self.run_command("central_symmetry_coloring")

    def change_aux_property_threshold(self, lower_upper, delta):
        if isinstance(lower_upper, int): lower_upper = str(lower_upper)
        self.run_command("change_aux_property_threshold %s %f" % (lower_upper, delta))

    def reset_aux_property_thresholds(self):
        self.run_command("reset_aux_property_thresholds")

    def toggle_aux_property_thresholds_saturation(self):
        self.run_command("toggle_aux_property_thresholds_saturation")

    def toggle_aux_property_thresholds_rigid(self):
        self.run_command("toggle_aux_property_thresholds_rigid")

    def rcut_patch(self, sym1, sym2, inc_dec, delta=None):
        self.run_command("rcut_patch start %s %s" % (sym1,sym2))
        if delta is None:
            self.run_command("rcut_patch %s" % inc_dec)
        else:
            self.run_command("rcut_patch %s %f" % (inc_dec, delta))
        self.run_command("rcut_patch finish")

    def select_gear(self, gear):
        self.run_command("select_gear %d" % gear)

    def cutting_plane(self, n, d, s):
        self.run_command("cutting_plane %d %f %f %f %f %f %f" % \
                                 (n, d[0], d[1], d[2], s[0], s[1], s[2]))

    def shift_cutting_plane_to_anchor(self, n):
        self.run_command("shift_cutting_plane_to_anchor %d" % n)

    def delete_cutting_plane(self, n):
        self.run_command("delete_cutting_plane %d" % n)

    def flip_cutting_plane(self, n):
        self.run_command("flip_cutting_plane %d" % n)

    def capture(self, filename, resolution=None):
        if resolution is None: resolution = ""
        format = filename[filename.rindex('.')+1:]
        self.run_command("capture %s %s %s" % (format, filename, resolution))

    def change_wireframe_mode(self, ):
        self.run_command("change_wireframe_mode")

    def change_cutting_plane_wireframe_mode(self):
        self.run_command("change_cutting_plane_wireframe_mode")

    def load_config(self, filename):
        self.run_command("load_config %s" % filename)

    def load_config_advance(self, command):
        self.run_command("load_config_advance %s" % command)

    def script_animate(self, filename):
        self.run_command("script_animate %s" % filename)

    def load_atom_color(self, filename):
        self.run_command("load_atom_color %s" % filename)

    def load_aux(self, filename):
        self.run_command("load_aux %s" % filename)

    def look_at_the_anchor(self):
        self.run_command("look_at_the_anchor")

    def observer_goto(self):
        self.run_command("observer_goto")

    def xtal_origin_goto(self, s):
        self.run_command("xtal_origin_goto %f %f %f" % (s[0], s[1], s[2]))

    def find_atom(self, i):
        if self.fortran_indexing: i = i-1
        self.run_command("find_atom %d" % i)

    def resize(self, width, height):
        self.run_command("resize %d %d" % (width, height))

    def change_aux_colormap(self, n):
        self.run_command("change_aux_colormap %d" % n)

    def print_atom_info(self, i):
        self.run_command("print_atom_info %d" % i)

    def save_atom_indices(self):
        self.run_command("save_atom_indices")

    def change_central_symm_neighbormax(self):
        self.run_command("change_central_symm_neighbormax")

    def timer(self, label):
        self.run_command("timer %s" % label)

    def isoatomic_reference_imprint(self):
        self.run_command("isoatomic_reference_imprint")

    def toggle_shell_viewer_mode(self):
        self.run_command("toggle_shell_viewer_mode")

    def toggle_xtal_mode(self):
        self.run_command("toggle_xtal_mode")

    def change_shear_strain_subtract_mean(self):
        self.run_command("change_shear_strain_subtract_mean")

    def draw_arrows(self, property, scale_factor=0.0, head_height=0.1, head_width=0.05, up=(0.0,1.0,0.0)):
        if property == 'off':
            self.run_command('draw_arrows off')
        else:
            self.run_command('draw_arrows %s %f %f %f %f %f %f' %
                             (str(property), scale_factor, head_height, head_width, up[0], up[1], up[2]))

    def wait(self):
        """Sleep until this AtomEye viewer has finished processing all queued events."""
        if not self.is_alive: 
            raise RuntimeError('is_alive is False')
        self._atomeye.wait(self._window_id)

views = {}
view = None

def show(obj, property=None, frame=0, window_id=None, nowindow=False, arrows=None, verbose=True, *arrowargs, **arrowkwargs):
    """Convenience function to show obj in the default AtomEye view

    If window_id is not None, then this window will be used. Otherwise
    the default window is used, initialising it if necessary.
    
    Returns instance of AtomEyeView."""
    global view

    # if window_id was passed in, then use that window
    if window_id is not None:
        views[window_id].show(obj, property, frame)
        return views[window_id]

    # otherwise use the default viewer, initialising it if necessary
    if view is None:
        if views.keys():
            view = views[views.keys()[0]]
            view.show(obj, property, frame)
        else:
            view = AtomEyeView(obj, property=property, frame=frame, nowindow=nowindow, arrows=arrows, verbose=verbose, *arrowargs, **arrowkwargs)
    else:
        view.show(obj, property, frame, arrows=arrows, *arrowargs, **arrowkwargs)

    return view
