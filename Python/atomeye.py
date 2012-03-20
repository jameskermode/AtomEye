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
import weakref
import atexit
from math import ceil, log10

__all__ = ['AtomEyeViewer', 'show']

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

    # remove directory at exit
    atexit.register(shutil.rmtree, path)
    
    return name, path


class MultipleExtMod(object):
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
        self._pkgname, self._pkgdir = _tmp_pkg()

        # add parent directory to sys.path if necessary
        if os.path.dirname(self._pkgdir) not in sys.path:
            sys.path.append(os.path.dirname(self._pkgdir))
        
        # copy the original module to the new package
        shutil.copy(srcpath, self._pkgdir)
        # import the module
        # __import__ returns the package, not the sub-module
        self._pkg = __import__(self._pkgname, globals(), locals(), [self.name])
        # return the module object
        self._module = getattr(self._pkg, self.name)
        # now add the module's stuff to this class instance
        self.__dict__.update(self._module.__dict__)

default_settings = {'n->xtal_mode': 1,
                    'n->suppress_printout': 1,
                    'n->bond_mode': 1,
                    'n->atom_r_ratio': 0.5,
                    'key->BackSpace': 'load_config_backward'
                    }

name_map = {'positions': 'pos',
            'masses'   : 'mass',
            'numbers'  : 'Z' }

viewers = {}

CONFIG_MAX_AUXILIARY = 64

class AtomEyeViewer(object):

    n_ext_modules = 0
    
    def __init__(self, atoms=None, viewer_id=None, copy=None, frame=0, delta=1,
                 nowindow=False, echo=False, block=False, verbose=True, **showargs):
        self.atoms = atoms
        self.frame = frame
        self.delta = delta
        self.echo = echo
        self.block = block
        self.fortran_indexing = False
        self.verbose = verbose
        self.is_alive = False
        self.selection = []

        if viewer_id is None:
            self.start(copy, nowindow)
        else:
            self._viewer_id = viewer_id
            self._window_id = viewer_id[1]
            self.is_alive = True
            viewers[self._viewer_id] = self

        self.show(**showargs)

    def _convert_atoms(self):
        sys.stderr.write('in convert_atoms.\n')
        self.current_atoms = None
        if self.atoms is None:
            title = '(null)'
            n_atom = 0
            cell = None
            arrays = None
        else:
            if hasattr(self.atoms, '__iter__'):
                sys.stderr.write('getting frame %d\n' % self.frame)
                self.current_atoms = self.atoms[self.frame]
                sys.stderr.write('done getting frame\n')
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

            if self.verbose:
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


    def _check_property_columns(self, auxprop):
        """
        Ensure auxprop is one of the first AUX_PROPERTY_COLORING
        columns in self.current_atoms by swapping it with an
        inessential property if necessary.
        """
            
        if not hasattr(self.current_atoms, 'properties'):
            return
        
        prop_idx = 0
        for key, value in self.current_atoms.properties.iteritems():
            ncols = len(value.shape) == 2 and value.shape[0] or 1
            prop_idx += ncols
            if key.lower() == auxprop.lower():
                break
        else:
            raise ValueError('Unknown Atoms property %s' % auxprop)

        if prop_idx >= CONFIG_MAX_AUXILIARY:
            for swapprop in self.current_atoms.properties:
                if swapprop.lower() not in ['pos', 'z', 'species']:
                    break
            self.current_atoms.properties.swap(auxprop, swapprop)

            
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
        self._atomeye = MultipleExtMod('_atomeye')
        self._atomeye.set_handlers(AtomEyeViewer.on_click,
                                   AtomEyeViewer.on_close,
                                   AtomEyeViewer.on_advance,
                                   AtomEyeViewer.on_new_window)

        self.is_alive = False
        self._window_id = self._atomeye.open_window(AtomEyeViewer.n_ext_modules,icopy,n_atom,cell,arrays,nowindow)
        self._viewer_id = (AtomEyeViewer.n_ext_modules, self._window_id)
        viewers[self._viewer_id] = self
        if AtomEyeViewer.n_ext_modules == 0:
            viewers['default'] = viewers[self._viewer_id]
        AtomEyeViewer.n_ext_modules += 1        
        while not self.is_alive:
            time.sleep(0.1)
        time.sleep(0.3)
        self._atomeye.set_title(self._window_id, title)
        self.update(default_settings)

    @staticmethod
    def on_click(mod, iw, idx):
        if (mod,iw) not in viewers:
            raise RuntimeError('Unexpected module id %d or window id %d' % (mod, iw))
        self = viewers[(mod,iw)]
        
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

    @staticmethod
    def on_advance(mod, iw, mode):
        if (mod,iw) not in viewers:
            raise RuntimeError('Unexpected window id %d' % iw)
        self = viewers[(mod, iw)]
        
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

    @staticmethod
    def on_close(mod, iw):
        print 'closing atomeye viewer', mod, iw
        if (mod, iw) not in viewers:
            raise RuntimeError('Unexpected window id %d' % iw)
        self = viewers[(mod,iw)]
        self.is_alive = False
        del viewers[self._viewer_id]
        if viewers['default'] is self:
            # check if we were the default viewer
            del viewers['default']
            if len(viewers) > 0:
                viewers['default'] = viewers.values()[0]

    @staticmethod
    def on_new_window(mod, iw):
        if (mod,iw) in viewers:
            viewers[(mod,iw)].is_alive = True
        else:
            new_viewer = AtomEyeViewer(viewer_id=(mod,iw))
            

    def show(self, obj=None, property=None, frame=None, arrows=None):
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

        sys.stderr.write('entering convert_atoms\n')
        title, n_atom, cell, arrays = self._convert_atoms() # also sets self.current_atoms
        sys.stderr.write('convert_atoms done\n')
                
        if property is not None and not isinstance(property,str) and hasattr(self.current_atoms, 'properties'):
            if isinstance(property,int):
                _show = [i == property for i in self.current_atoms.indices]
            elif isinstance(property, list) or isinstance(property, tuple) or isinstance(property, set):
                _show = [i in property for i in self.current_atoms.indices]
            else:
                _show = property

            if self.current_atoms.has_property('_show'):
                self.current_atoms.remove_property('_show')
            self.current_atoms.add_property('_show', _show)
            property = '_show'

        if self.current_atoms is not None:
            sys.stderr.write('calling load_atoms with current_atoms=%r\n' % self.current_atoms)

            self._atomeye.load_atoms(self._window_id, title, n_atom, cell, arrays)
            if property is not None:
                self.aux_property_coloring(property)
            if arrows is not None:
                self.draw_arrows(arrows)

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
        self._check_property_columns(auxprop)
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

    def first(self):
        self.run_command("load_config_first")

    def last(self):
        self.run_command("load_config_last")

    def forward(self):
        self.run_command("load_config_forward")

    def backward(self):
        self.run_command("load_config_backward")

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
        if property is None:
            self.run_command('draw_arrows off')
        else:
            self._check_property_columns(property)
            self.run_command('draw_arrows %s %f %f %f %f %f %f' %
                             (str(property), scale_factor, head_height, head_width, up[0], up[1], up[2]))

    def wait(self):
        """Sleep until this AtomEye viewer has finished processing all queued events."""
        if not self.is_alive: 
            raise RuntimeError('is_alive is False')
        self._atomeye.wait(self._window_id)

    def get_visible(self):
        """Return list of indices of atoms currently visible in this viewer."""
        indices = self._atomeye.get_visible()
        if numpy.any(indices > len(self.current_atoms)):
            # atoms duplicated due to narrow cell (n->small_cell_err_handler == 1)
            indices = list(set([idx % len(self.current_atoms) for idx in indices ]))
        if self.fortran_indexing:
            indices = [idx+1 for idx in indices]
        return indices


def show(obj, property=None, frame=0, viewer=None,
         nowindow=False, newwindow=False, arrows=None, verbose=True):
    """Show `obj` with AtomEye, opening a new window if necessary.

    The viewer to be used is the first item on the following list which is defined:

      1. `viewer` argument - should be an AtomEyeViewer instance
      3. `obj.viewer` attribute - should be a weak reference to an AtomEyeViewer instance
      4. The default viewer - `atomeye.viewers["default"]`

    If no viewers exist, or if newwindow=True, a new viewer is created.
    In all cases, the instance of AtomEyeViewer used is returned."""

    if not newwindow:
        # if obj has been viewed before, viewer will have been saved with weak reference
        if viewer is None and hasattr(obj, 'viewer') and type(obj.viewer) is weakref.ref:
            viewer = obj.viewer() 

        if viewer is None and 'default' in viewers:
            # Use default (i.e. first created) viewer
            viewer = viewers['default']

    if viewer is None:
        # We need to create a new viewer
        viewer = AtomEyeViewer(obj,
                               verbose=verbose,
                               nowindow=nowindow,
                               property=property,
                               frame=frame,
                               arrows=arrows)
    else:
        viewer.show(obj,
                    property=property,
                    frame=frame,
                    arrows=arrows)

    obj.viewer = weakref.ref(viewer) # save a reference to viewer for next time
    return viewer


