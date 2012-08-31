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
import atexit
from math import ceil, log10

__all__ = ['AtomEyeViewer']

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

    count = {}

    def __init__(self, name):
        self.count = MultipleExtMod.count[name] = MultipleExtMod.count.setdefault(name, 0) + 1
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


default_state = {
    'variables' : {'n->xtal_mode': 1,
                   'n->suppress_printout': 1,
                   'n->bond_mode': 1,
                   'n->atom_r_ratio': 0.5,
                   'key->BackSpace': 'load_config_backward'
                   },
    'commands': ['xtal_origin_goto 0.5 0.5 0.5',
                 'toggle_parallel_projection'],
    'rcut_patches': [('Si', 'Si', 0.5)]
}

name_map = {'positions': 'pos',
            'masses'   : 'mass',
            'numbers'  : 'Z' }

viewers = {}

class AtomEyeViewer(object):
    """
    View an atomic configuration or trajectory with AtomEye.
    """

    CONFIG_MAX_AUXILIARY = 64
    
    def __init__(self, atoms=None, viewer_id=None, copy=None, frame=0, delta=1,
                 nowindow=False, echo=False, block=False, verbose=True, fortran_indexing=False,
                 **showargs):
        """
        Construct a new AtomEye viewer window. Each viewer class communicates with one
        AtomEye thread.
        
        `atoms` - configuration or trajectory to view. Should be a :class:`quippy.Atoms`
                  or :class:`ase.Atoms` instance, or a list of such instances.
        `viewer_id` - if None, open a new viewer. Otherwise call the :meth:`show`() method in the existing viewer with this ID.
        `copy` - viewer ID of another viewer from which to copy the viewpoint and other default settings
        `frame` - initial frame to show (should be in range 0..len(atoms)-1)
        `delta` - increment/decement rate for frames when [Insert] and [Delete] are pressed
        `nowindow` - if True, open AtomEye without a visible window. Useful for faster rendering of movies
        `echo` - if True, echo all commands to the screen. Good for debugging.
        `block` - if True, wait for commands to finish executing in AtomEye thread before returning (i.e. run asynchronously)
        `verbose` - if True, print information when changing frame and when an atom is clicked

        Additional keyword arguments are passed along to the :method:`show`() method.
        """
        self.atoms = atoms
        self._current_atoms = None
        self._previous_atoms = None
        self.frame = frame
        self.delta = delta
        self.echo = echo
        self.block = block
        self.fortran_indexing = fortran_indexing
        self.verbose = verbose
        self.is_alive = False

        if viewer_id is None:
            self.start(copy, nowindow)
        else:
            self._viewer_id = viewer_id
            self._window_id = viewer_id[1]
            self.is_alive = True
            viewers[self._viewer_id] = self

        self.show(**showargs)

           
    def start(self, copy=None, nowindow=False):
        """
        Start the AtomEye thread, wait for it to load and apply default commands.
        """
        if self.is_alive: return
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
                                   AtomEyeViewer.on_new_window,
                                   AtomEyeViewer.on_redraw)

        self.is_alive = False
        self._module_id = self._atomeye.count
        self._window_id = len([ viewer for viewer in viewers if viewer[0] == self._module_id ])
        self._viewer_id = (self._module_id, self._window_id)
        print 'Initialising AtomEyeViewer with module_id %d and window id %s' % self._viewer_id
        viewers[self._viewer_id] = self
        atomeye_window_id = self._atomeye.open_window(self._module_id, icopy, nowindow)
        #assert atomeye_window_id == self._window_id
        while not self.is_alive:
            time.sleep(0.1)
        time.sleep(0.3)
        self.set_state(default_state)
        self.wait()

    def _click_hook(self, at, idx):
        pass

    def _enter_hook(self, at):
        pass

    def _exit_hook(self, at):
        pass

    def _close_hook(self):
        pass

    def _redraw_hook(self, at):
        pass

    def _property_hook(self, at, auxprop):
        return auxprop

    @staticmethod
    def on_click(mod, iw, idx):
        if (mod,iw) not in viewers:
            raise RuntimeError('Unexpected module id %d or window id %d' % (mod, iw))
        self = viewers[(mod,iw)]
        at = self.gca()
        
        if at is None:
            return
        if idx >= len(at):
            idx = idx % len(at)
        if self.fortran_indexing:
            idx = idx + 1 # atomeye uses zero based indices
        self._click_hook(at, idx)

    @staticmethod
    def on_advance(mod, iw, mode):
        if (mod,iw) not in viewers:
            raise RuntimeError('Unexpected window id %d' % iw)
        self = viewers[(mod, iw)]
        if not hasattr(self.atoms,'__iter__'): return
        if mode not in ['forward', 'backward', 'first', 'last']:
            raise RuntimeError('Unexpected advance mode "%s"' % mode)
        getattr(self, mode)()

    @staticmethod
    def on_close(mod, iw):
        if (mod, iw) not in viewers:
            raise RuntimeError('Unexpected window id %d' % iw)
        self = viewers[(mod,iw)]
        self.is_alive = False
        del viewers[self._viewer_id]
        self._close_hook()

    @staticmethod
    def on_new_window(mod, iw):
        if (mod,iw) in viewers:
            viewers[(mod,iw)].is_alive = True
        else:
            new_viewer = AtomEyeViewer(viewer_id=(mod,iw))
            
    @staticmethod
    def on_redraw(mod, iw):
        if (mod, iw) not in viewers:
            raise RuntimeError('Unexpected window id %d' % iw)
        self = viewers[(mod,iw)]

        # keep a reference to old atoms around so memory doesn't get free'd prematurely
        self._previous_atoms = self._current_atoms
        if self._previous_atoms is not None:
            self._exit_hook(self._previous_atoms)
        self._current_atoms = None
        
        if self.atoms is None:
            title = '(null)'
            n_atom = 0
            cell = None
            arrays = None
            redraw = 0
        else:
            name = 'Atoms'
            if hasattr(self.atoms, 'filename'):
                name = self.atoms.filename
            if hasattr(self.atoms, 'name'):
                name = self.atoms.name
            self._current_atoms = self.gca(update=True)
            if hasattr(self.atoms, '__iter__'):
                fmt = "%%0%dd" % ceil(log10(len(self.atoms)+1))
                title = '%s frame %s length %s' % (name, fmt % self.frame, fmt % len(self.atoms))
            else:
                title = name

            self._enter_hook(self._current_atoms)
            n_atom = len(self._current_atoms)
            try:
                self.fortran_indexing = self._current_atoms.fortran_indexing
            except AttributeError:
                self.fortran_indexing = False

            cell = self._current_atoms.get_cell()
            pbc = self._current_atoms.get_pbc()
            pos = self._current_atoms.positions
            
            for i, p in enumerate(pbc):
                if not p:
                    cell[i,i] = max(1.0, 2*(pos[:,i].max()-pos[:,i].min()))

            try:
                arrays = self._current_atoms.properties
            except AttributeError:
                arrays = {}
                for key,value in self._current_atoms.arrays.iteritems():
                    arrays[name_map.get(key,key)] = value

            redraw = 1 # FIXME, we should decide if we really have to redraw here

        if redraw and self.verbose:
            print '-'*len(title)
            print title
            print '-'*len(title)
            print 'Number of atoms: %d' % n_atom
            print 'Fortran indexing: %r' % self.fortran_indexing
            print 'Unit cell:'
            print cell
            self._redraw_hook(self._current_atoms)
            print '\n'
            sys.stdout.flush()
        
        return (redraw, title, n_atom, cell, arrays)

    def gcat(self, update=False):
        """Get current atoms - return Atoms object currently being visualised.

        If update=False (the default), we return what is currently being visualised,
        even if this is not in sync with self.atoms[self.frame]."""

        if not update and self._current_atoms is not None:
            return self._current_atoms
        else:
            if hasattr(self.atoms, '__iter__'):
                return self.atoms[self.frame % len(self.atoms)]
            else:
                return self.atoms

    def scat(self, atoms, frame=None):
        """Set current atoms (and optionally also current frame)"""
        if atoms is not None:
            self.atoms = atoms
        if frame is not None and hasattr(self.atoms, '__iter__'):
            self.frame = frame % len(self.atoms)

    def show(self, atoms=None, property=None, frame=None, arrows=None):
        """
        Update what is shown in this AtomEye viewer window.

        `atoms` should be a quippy.Atoms or ase.Atoms instance, or a list of instances.
        `property` should be the name of the auxiliary property used to colour the atoms (e.g. "charge")
        `frame` is the (zero-based) index of the frame to show.
        `arrows` is the name of a vector property to use to draw arrows on the atoms (e.g. "force")

        When called with no arguments, show() is equivalent to redraw()
        """
        
        if not self.is_alive:
            raise RuntimeError('is_alive is False')
        self.sca(atoms, frame)
        if property is not None:
            self.aux_property_coloring(property)
        if arrows is not None:
            self.draw_arrows(arrows)
        self.redraw()
                
    def redraw(self):
        """
        Redraw this AtomEye window, keeping Atoms and settings the same
        """
        self._atomeye.redraw(self._window_id)

    def run_command(self, command):
        """
        Run a command in this AtomEye thread. The command is queued for later execution.
        """
        if not self.is_alive: 
            raise RuntimeError('is_alive is False')
        if self.echo:
            print command.strip()
        try:
            self._atomeye.run_command(self._window_id, command)
        except RuntimeError as err:
            if str(err) == 'atomeyelib_queueevent: too many atomeyelib events.':
                self.wait()
                self._atomeye.run_command(self._window_id, command)
            else:
                raise
        if self.block:
            self.wait()

    def run_script(self, script):
        """
        Run commands from the file script, in a blocking fashion.
        """
        if type(script) == type(''):
            script = open(script)
            
        for line in script:
            self.run_command(line)
            self.wait()

    def __call__(self, command):
        self.run_command(command)

    def close(self):
        """
        Close this viewer window.
        """
        self.run_command('close')

    def setp(self, key, value):
        """
        Run the AtomEye command "set key value".
        """
        self.run_command("set %s %s" % (str(key), str(value)))

    def save(self, filename):
        """
        Save AtomEye viewer settings to a file.
        """
        self.run_command("save %s" % str(filename))

    def update(self, state):
        """
        Update settings from the dictionary D.
        """
        for k, v in state.iteritems():
            self.setp(k, v)

    def set_state(self, state):
        for key, value in state.iteritems():
            if key == 'variables':
                self.update(value)
            elif key == 'commands':
                for command in value:
                    self.run_command(command)
            elif key == 'rcut_patches':
                for (sym1, sym2, rcut) in value:
                    self.rcut_patch(sym1, sym2, float(rcut))
            else:
                setattr(self, key, value)
        self.redraw()
                
    def get_state(self):
        fd, name = tempfile.mkstemp()
        os.close(fd)
        self.save(name)
        self.wait()
        time.sleep(0.2)
        fh = open(name, 'r')
        lines = fh.readlines()
        variables = {}
        commands = []
        for line in lines:
            line = line.strip()
            if line.startswith('set'):
                dummy, key, value = line.split(None,2)
                variables[key] = value
            else:
                commands.append(line)
        fh.close()
        os.unlink(name)
        state = {}
        state['variables'] = variables
        state['commands'] = commands

        odict = self.__dict__.copy()
        for key in odict.keys():
            if key.startswith('_'):
                del odict[key]
        del odict['atoms']
        state.update(odict)

        return state

    def __getstate__(self):
        return self.get_state()

    def __setstate__(self, state):
        self.set_state(state)

    def load_script(self, filename):
        """
        Load AtomEye viewer settings from a file - run_script() is more robust as it's blocking.
        """
        self.run_command("load_script %s" % str(filename))

    def key(self, key):
        """
        Simulate pressing `key` on the keyboard.
        """
        self.run_command("key %s" % key)

    def toggle_coordination_coloring(self):
        """
        Turn on or off colouring by coordination number (key "k")
        """
        self.run_command("toggle_coordination_coloring")

    def translate(self, axis, delta):
        """
        Translate system along `axis` by an amount `delta` (key "Ctrl+left/right/up/down")
        """
        self.run_command("translate %d %f " % (axis, delta))

    def shift_xtal(self, axis, delta):
        """
        Shift crystal within periodic boundaries along `axis` by `delta` (key "Shift+left/right/up/down").
        """
        self.run_command("shift_xtal %d %f" % (axis, delta))

    def rotate(self, axis, theta):
        """
        Rotate around `axis` by angle `theta`.
        """
        self.run_command("rotate %d %f" % (axis, theta))

    def advance(self, delta):
        """
        Move the camera forward by `delta`.
        """
        self.run_command("advance %f" % delta)

    def shift_cutting_plane(self, delta):
        """
        Shift the current cutting plane by an amount `delta`.
        """
        self.run_command("shift_cutting_plane %f" % delta)

    def change_bgcolor(self, color):
        """
        Change the viewer background colour to `color`, which should be a RGB tuple with three floats in range 0..1.
        """
        self.run_command("change_bgcolor %f %f %f" % (color[0], color[1], color[2]))

    def change_atom_r_ratio(self, delta):
        """
        Change the size of the balls used to draw the atoms by `delta`.
        """
        self.run_command("change_atom_r_ratio %f" % delta)

    def change_bond_radius(self, delta):
        """
        Change the radius of the cylinders used the draw bonds by `delta`.
        """
        self.run_command("change_bond_radius %f" % delta)

    def change_view_angle_amplification(self, delta):
        self.run_command("change_view_angle_amplification %f" % delta)

    def toggle_parallel_projection(self):
        """
        Toggle between parallel and perspective projections.
        """
        self.run_command("toggle_parallel_projection")

    def toggle_bond_mode(self):
        """
        Turn on or off bonds.
        """
        self.run_command("toggle_bond_mode" )

    def toggle_small_cell_mode(self):
        """
        Toggle between two different behaviours for when cell is smaller than r_cut/2:
         1. clip cell - some neigbours may be lost (default)
         2. replicate cell along narrow directions
        """
        self.run_command("toggle_small_cell_mode")
        self.redraw()

    def normal_coloring(self):
        """
        Return to normal colouring of the atoms (key "o").
        """
        self.run_command("normal_coloring")

    def aux_property_coloring(self, auxprop):
        """
        Colour the atoms by the auxiliary property with name or index `auxprop`.
        """
        auxprop = self._property_hook(self.gca(), auxprop)
        self.run_command("aux_property_coloring %s" % str(auxprop))

    def central_symmetry_coloring(self):
        """
        Colour atoms by centro-symmetry parameter.
        """
        self.run_command("central_symmetry_coloring")

    def change_aux_property_threshold(self, lower, upper):
        """
        Change the lower and upper aux property thresholds.
        """
        self.run_command("change_aux_property_threshold lower %f" % lower)
        self.run_command("change_aux_property_threshold upper %f" % upper)

    def reset_aux_property_thresholds(self):
        """
        Reset aux property thresholds to automatic values.
        """
        self.run_command("reset_aux_property_thresholds")

    def toggle_aux_property_thresholds_saturation(self):
        """
        Toggle between saturated colouring and invisibility for values outside aux prop thresholds.
        """
        self.run_command("toggle_aux_property_thresholds_saturation")

    def toggle_aux_property_thresholds_rigid(self):
        """
        Toggle between floating and rigid aux property thresholds when moving between frames
        """
        self.run_command("toggle_aux_property_thresholds_rigid")

    def rcut_patch(self, sym1, sym2, delta):
        """
        Change the cutoff distance for `sym1`--`sym2` bonds by `delta.

        e.g. to increase cutoff for Si-Si bonds by 0.5 A use
             viewer.rcut_patch('Si', 'Si', 0.5)
        """
        self.run_command("rcut_patch start %s %s" % (sym1,sym2))
        self.run_command("rcut_patch %s" % str(delta))
        self.run_command("rcut_patch finish")

    def select_gear(self, gear):
        self.run_command("select_gear %d" % gear)

    def cutting_plane(self, n, d, s):
        """
        Create a new cutting plane with index `n`, normal `d`, and fractional displacement `s`.
        """
        self.run_command("cutting_plane %d %f %f %f %f %f %f" % \
                                 (n, d[0], d[1], d[2], s[0], s[1], s[2]))

    def shift_cutting_plane_to_anchor(self, n):
        self.run_command("shift_cutting_plane_to_anchor %d" % n)

    def delete_cutting_plane(self, n):
        self.run_command("delete_cutting_plane %d" % n)

    def flip_cutting_plane(self, n):
        self.run_command("flip_cutting_plane %d" % n)

    def capture(self, filename, resolution=None):
        """
        Save current view to `filename`. Format is determined from file extension: .png, .jpeg, or .eps.
        """
        if resolution is None: resolution = ""
        format = filename[filename.rindex('.')+1:]
        self.run_command("capture %s %s %s" % (format, filename, resolution))

    def change_wireframe_mode(self, ):
        self.run_command("change_wireframe_mode")

    def change_cutting_plane_wireframe_mode(self):
        self.run_command("change_cutting_plane_wireframe_mode")

    def frame(self, frame):
        """
        Set the current frame to `frame`.
        """
        
        self.frame = frame % len(self.atoms)
        self.redraw()

    def first(self):
        """
        Show the fisrt frame (frame 0).
        """

        self.frame = 0
        self.redraw()

    def last(self):
        """
        Show the last frame, i.e. len(self.atoms)-1
        """

        self.frame = len(self.atoms)-1
        self.redraw()

    def forward(self, delta=None):
        """
        Move forward by `delta` frames (default value is self.delta).
        """

        delta = delta or self.delta
        self.frame += delta
        self.frame = self.frame % len(self.atoms)
        self.redraw()

    def backward(self, delta=None):
        """
        Move backward by `delta` frames (default values is self.delta).
        """

        delta = delta or self.delta
        self.frame -= delta
        self.frame = self.frame % len(self.atoms)
        self.redraw()

    def load_atom_color(self, filename):
        """
        Load atom colours from a .clr file.
        """
        self.run_command("load_atom_color %s" % filename)

    def load_aux(self, filename):
        """
        Load aux property values from a .aux file.
        """
        self.run_command("load_aux %s" % filename)

    def look_at_the_anchor(self):
        self.run_command("look_at_the_anchor")

    def observer_goto(self):
        self.run_command("observer_goto")

    def xtal_origin_goto(self, s):
        """
        Move the crystal origin to fractional coordinates s (e.g. [0.5, 0.5, 0.5).
        """
        self.run_command("xtal_origin_goto %f %f %f" % (s[0], s[1], s[2]))

    def find_atom(self, i):
        """
        Set the anchor to the atom with index `i`.
        """
        if self.fortran_indexing: i = i-1
        self.run_command("find_atom %d" % i)

    def resize(self, width, height):
        """
        Resize the current window to `width` x `height` pixels.
        """
        
        self.run_command("resize %d %d" % (width, height))

    def change_aux_colormap(self, n):
        """
        Select the `n`-th auxiliary property colourmap. 
        """
        self.run_command("change_aux_colormap %d" % n)

    def draw_arrows(self, property, scale_factor=0.0, head_height=0.1, head_width=0.05, up=(0.0,1.0,0.0)):
        """
        Draw arrows on each atom, based on the vector property with name `property` (e.g. "force").

        Passing `property`=None turns off previous arrows. `scale_factor` can be used to change length of
        arrows (1 unit = 1 Angstrom; default value of 0.0 means autoscale). `head_height` and  `head_width`
        set height and with of arrow heads in Angstrom. The `up` vector controls the plane in which arrow heads
        are drawn (default [0.,1.,0.]).
        """
        
        if property is None:
            self.run_command('draw_arrows off')
        else:
            property = self._property_hook(self.gca(), property)
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
        at = self.gca()
        if numpy.any(indices > len(at)):
            # atoms duplicated due to narrow cell (n->small_cell_err_handler == 1)
            indices = list(set([idx % len(at) for idx in indices ]))
        if self.fortran_indexing:
            indices = [idx+1 for idx in indices]
        return indices


