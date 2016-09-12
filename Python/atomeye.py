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

"""

This module provides the :class:`AtomEyeViewer` class, which is a high-level
interface for interactive visualisation of Atoms objects using a `modified version
<http://www.jrkermode.co.uk/AtomEye>`_ of the `AtomEye <http://mt.seas.upenn.edu/Archive/Graphics/A/>`_
atomistic configuration viewer.

:class:`~quippy.atoms.Atoms` and :class:`~quippy.io.AtomsList`
objects can also be visualised with the :mod:`qlab` module.

For example, to create and visualise and 8-atom silicon bulk cell::

   from quippy.structures import diamond, supercell
   from atomeye import AtomEyeViewer

   d = diamond(5.43, 14)
   viewer = AtomEyeViewer(d)

A window will pop up containing the silicon unit cell, which you can manipulate
with the mouse, by issuing commands in the the python console or
with a variety of AtomEye `keyboard shortcuts
<http://li.mit.edu/Archive/Graphics/A3/A3.html#keybinding>`_. To save an image in PNG format::

   viewer.capture('si8.png')

.. image:: si8.png
   :align: center

Then, to change the display to a :math:`2\times2\times2` supercell of
bulk silicon, change the background colour to black, set the size and
save an image you'd do the following::

   at = supercell(d, 2, 2, 2)
   viewer.show(at)

   viewer.change_bgcolor((0, 0, 0))
   viewer.resize(400,300)
   viewer.capture('si2x2x2.png')

.. image:: si2x2x2.png
   :align: center

Module attributes:

.. attribute:: viewers

   Dictionary mapping `window_id` to :class:`AtomEyeView` instances. There
   is one entry for each currently open AtomEye view window.

.. attribute:: default_state

   Dictionary of key/value pairs passed to :meth:`AtomEyeView.update` when
   a new window is created. Change this dictionary to modfiy the properties
   of new windows. The initial value is as follows::

     default_state = {
	 'variables' : {'n->xtal_mode': 1,
			'n->suppress_printout': 1,
			'n->bond_mode': 1,
			'n->atom_r_ratio': 0.5,
			'n->small_cell_err_handler': 1,
			'key->BackSpace': 'load_config_backward'
			},
	 'commands': ['xtal_origin_goto 0.5 0.5 0.5',
		      'toggle_parallel_projection'],
	 'cutoff_lengths': []
     }
"""

import sys
import time
import imp
import os
import os.path
import shutil
import sys
import tempfile
import atexit
from math import ceil, log10
import numpy as np

__all__ = ['AtomEyeViewer', 'view', 'redraw', 'run_command', 'run_script', 'close',
           'setp', 'save_script', 'toggle_coordination_coloring', 'translate',
           'shift_xtal', 'rotate', 'advance', 'shift_cutting_plane', 'change_bgcolor',
           'change_atom_r_ratio', 'change_bond_radius', 'change_view_angle_amplification',
           'toggle_parallel_projection', 'toggle_bond_mode', 'toggle_small_cell_mode',
           'normal_coloring', 'aux_property_coloring', 'central_symmetry_coloring',
           'change_aux_property_threshold', 'reset_aux_property_thresholds',
           'toggle_aux_property_thresholds_saturation', 'toggle_aux_property_thresholds_rigid',
           'rcut_patch', 'select_gear', 'cutting_plane', 'shift_cutting_plane_to_anchor',
           'delete_cutting_plane', 'flip_cutting_plane', 'capture', 'change_wireframe_mode',
           'change_cutting_plane_wireframe_mode', 'get_frame', 'set_frame', 'get_delta',
           'set_delta', 'first', 'last', 'forward', 'backward', 'load_atom_color',
           'load_aux', 'look_at_the_anchor', 'observer_goto', 'xtal_origin_goto',
           'find_atom', 'resize', 'change_aux_colormap', 'draw_arrows', 'wait', 'display']

try:
    # try to use same fortran_indexing setting as quippy, if it's installed
    from quippy import get_fortran_indexing, set_fortran_indexing
    
except ImportError:
    # quippy is not available, so we define our own fortran_indexing setting

    _fortran_indexing = False
    
    def get_fortran_indexing():
        global _fortran_indexing
        return _fortran_indexing

    def set_fortran_indexing(fortran_indexing):
        global _fortran_indexing
        _fortran_indexing = fortran_indexing

    __all__.extend(['get_fortran_indexing', 'set_fortran_indexing'])
                             

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
                   'n->small_cell_err_handler': 1,
                   'key->BackSpace': 'load_config_backward'
                   },
    'commands': ['xtal_origin_goto 0.5 0.5 0.5',
                 #'toggle_parallel_projection',
                 'resize 800 600',
                 'change_atom_r_ratio -2.0'],
    'cutoff_lengths': []
}

name_map = {'positions': 'pos',
            'masses'   : 'mass',
            'numbers'  : 'Z' }

viewers = {}

class AtomEyeViewer(object):
    """
    View an atomic configuration or trajectory with AtomEye

    Class to represent an AtomEye viewer window. Each viewer class
    communicates with one AtomEye thread.
 
    There are wrapper methods for most of the `AtomEye 3 commands
    <http://mt.seas.upenn.edu/Archive/Graphics/A3/A3.html#commands>`_.
    The names of these methods match the names of the correspondning
    commands, and the arguments follow the syntax given on the AtomEye
    3 help page.
    
    Additional keyword arguments are passed along to the
    :meth:`show` method.

    Import :class:`AtomEyeView` attributes:

    .. attribute:: atoms
 
       :class:`Atoms` object or sequence being viewed. This will be set
       to ``None`` if this instance was created without an ``obj``
       parameter, which means we're viewing the ``A3`` logo.
 
    .. attribute:: frame
 
       Current frame, in range 1 to `len(self.atoms)`.
 
    .. attribute:: delta
 
       Frame increment rate when :kbd:`Delete` and :kbd:`Insert` are
       preseed. Equivalent to AtomEye ``n->glob_advance`` setting.
 
    .. attribute:: echo
 
       If set to True, echo all AtomEye commands to stdout
 
    .. attribute:: block
 
       If set to True, wait for all AtomEye command to finish
       executing before returning from function calls. Default is
       False.
 
    .. attribute:: verbose
 
       If set to True (default), print frame paramters on each
       :meth:`redraw`, and print information about each atom when it is
       right clicked.

    Parameters
    ----------

    atoms     : :class:`quippy.atoms.Atoms` or :class:`ase.atoms.Atoms` object, or a list of objects.
       Configuration or trajectory to view.
    viewer_id : integer or None
       If None, open a new viewer. Otherwise call the
       :meth:show() method in the existing viewer with this ID.
    copy      : integer or None
       Viewer ID of another viewer from which to copy the viewpoint
       and other default settings.
    frame     : integer
       Initial frame to show (should be in range ``0..len(atoms)-1``)
    delta     : integer
       Increment/decrement rate for frames when [Insert] and [Delete]
       are pressed
    nowindow  : bool
       If True, open AtomEye without a visible window. Useful
       for faster rendering of movies
    echo      : bool
       If True, echo all commands to the screen.
       Useful mainly for debugging.
    block     : bool
       If True, wait for commands to finish executing in AtomEye
       thread before returning (i.e. run asynchronously)
    verbose   : bool
       If True, print information when changing frame and
       when an atom is clicked

    """

    CONFIG_MAX_AUXILIARY = 64
    
    def __init__(self, atoms=None, viewer_id=None, copy=None, frame=0, delta=1,
                 nowindow=False, echo=False, block=False, verbose=False,
                 **showargs):
        self.atoms = atoms
        self._current_atoms = None
        self._previous_atoms = None
        self._frame = frame
        self.delta = delta
        self.echo = echo
        self.block = block
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
        Start the AtomEye thread, wait for it to load and apply
        default commands.

        Parameters
        ----------

        copy     : integer
           Viewer ID of another AtomEye window from which to
           copy settings
        nowindow : bool
           If True, create an AtomEye thread without a visible
           window. Useful for rendering movies.
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
        at = self.gcat()
        
        if at is None:
            return
        if idx >= len(at):
            idx = idx % len(at)
        if get_fortran_indexing():
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
            if hasattr(self._current_atoms, 'filename') and self._current_atoms.filename is not None:
                name = self._current_atoms.filename
            if hasattr(self._current_atoms, 'name') and self._current_atoms.name is not None:
                name = self._current_atoms.name
            self._current_atoms = self.gcat(update=True)
            if hasattr(self.atoms, '__iter__'):
                fmt = "%%0%dd" % ceil(log10(len(self.atoms)+1))
                title = '%s frame %s length %s' % (name, fmt % self._frame, fmt % len(self.atoms))
            else:
                title = name

            self._enter_hook(self._current_atoms)
            n_atom = len(self._current_atoms)

            cell = self._current_atoms.get_cell()
            pbc = self._current_atoms.get_pbc()
            pos = self._current_atoms.positions
            
            for i, p in enumerate(pbc):
                if not p:
                    cell[i,i] = max(cell[i,i], max(1.0, 2*(pos[:,i].max()-pos[:,i].min())))

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
            print 'Fortran indexing: %r' % get_fortran_indexing()
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
                return self.atoms[self._frame % len(self.atoms)]
            else:
                return self.atoms

    def scat(self, atoms, frame=None):
        """Set current atoms (and optionally also current frame)"""
        if atoms is not None:
            self.atoms = atoms
        if frame is not None and hasattr(self.atoms, '__iter__'):
            self._frame = frame % len(self.atoms)
        self.redraw()

    def show(self, atoms=None, property=None, frame=None, arrows=None):
        """
        Update what is shown in this AtomEye viewer window.

        When called with no arguments, :meth:`show` is equivalent to
        :meth:`redraw`.

        Parameters
        ----------

        atoms    : class:`quippy.atoms.Atoms or `ase.atoms.Atoms`
                   instance, or a list of instances
        property : name of the quippy :``~quippy.atoms.properties``
                   entry or ASE :attr:`ase.atoms.arrays` entry used to
                   colour the atoms (e.g. ``"charge"``)
        frame    : Zero-based index of the frame to show
                   (applies only when `atoms` is a list of Atoms objects)
        arrows   : is the name of a vector property to use to draw arrows
                   on the atoms (e.g. ``"force"``)
        """
        
        if not self.is_alive:
            raise RuntimeError('is_alive is False')
        self.scat(atoms, frame)
        if property is not None:
            self.aux_property_coloring(property)
        if arrows is not None:
            self.draw_arrows(arrows)
        self.redraw()
                
    def redraw(self):
        """
        Redraw this AtomEye window, keeping Atoms and settings the same.
        """
        self._atomeye.redraw(self._window_id)

    def run_command(self, command):
        """
        Run a command in this AtomEye thread.

        The command is queued for later execution, unless :attr:`block` is True.

        This functionality is also available by calling an instance
        directly, i.e. the following commands are equivalent::
      
	   viewer.run_command('toggle_coordination_coloring')
           viewer('toggle_coordination_coloring')

        Parameters
        ----------

        command : string
           The command to pass to AtomEye
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
        Update settings from the dictionary `state`.

        Runs the AtomEye command ``set key value`` for each pair. Valid settings
        are listed on the `AtomEye 3 settings help page
        <http://mt.seas.upenn.edu/Archive/Graphics/A3/A3.html#redraw>`_
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
            elif key == 'cutoff_lengths':
                for (sym1, sym2, cutoff) in value:
                    self.rcut_patch(sym1, sym2, float(cutoff), absolute=True)
            else:
                setattr(self, key, value)
            self.wait()
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
        Load AtomEye viewer settings from a file using the AtomEye ``load_script`` command

        :meth:`run_script` is more robust as the script is run line by
        line in a blocking sense.
        """
        self.run_command("load_script %s" % str(filename))

    def key(self, key):
        """
        Simulate pressing `key` on the keyboard.

        The syntax for keystrokes is described on the `AtomEye 3
        commands help page
        <http://mt.seas.upenn.edu/Archive/Graphics/A3/A3.html#commands>`_
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
        r, g, b = color
        self.run_command("change_bgcolor %f %f %f" % (r, g, b))

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
        auxprop = self._property_hook(self.gcat(), auxprop)
        self.redraw() # ensure auxprop is available
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

    def rcut_patch(self, sym1, sym2, value, absolute=False):
        """
        Change the cutoff distance for `sym1`--`sym2` bonds by `delta`.

        e.g. to increase cutoff for Si-Si bonds by 0.5 A use::

             viewer.rcut_patch('Si', 'Si', 0.5)

        With `absolute` set to True, `value` is used to set the
        absolute cutoff distance for `sym1`--`sym2` bonds, e.g.::

             viewer.rcut_patch('Si', 'Si', 2.50, True)
        """
        self.run_command("rcut_patch start %s %s" % (sym1,sym2))
        cmd = "rcut_patch %s" % value
        if absolute:
            cmd += " absolute"
        self.run_command(cmd)
        self.run_command("rcut_patch finish")

    def select_gear(self, gear):
        self.run_command("select_gear %d" % gear)

    def cutting_plane(self, n, d, s):
        """
        Create a new cutting plane with index `n`, normal `d`, and fractional displacement `s`.
        """
        da, db, dc = d
        sa, sb, sc = s
        self.run_command("cutting_plane %d %f %f %f %f %f %f" % \
                                 (n, da, db, dc, sa, sb, sc))

    def shift_cutting_plane_to_anchor(self, n):
        self.run_command("shift_cutting_plane_to_anchor %d" % n)

    def delete_cutting_plane(self, n):
        self.run_command("delete_cutting_plane %d" % n)

    def flip_cutting_plane(self, n):
        self.run_command("flip_cutting_plane %d" % n)

    def capture(self, filename, resolution=None):
        """
        Render the current view to image `filename`

        Format is determined from file extension: .png, .jpeg, or .eps.
        """
        if resolution is None: resolution = ""
        format = filename[filename.rindex('.')+1:]
        self.wait()
        self.run_command("capture %s %s %s" % (format, filename, resolution))

    def change_wireframe_mode(self, ):
        self.run_command("change_wireframe_mode")

    def change_cutting_plane_wireframe_mode(self):
        self.run_command("change_cutting_plane_wireframe_mode")

    def get_frame(self):
        return self._frame

    def set_frame(self, frame):
        self._frame = frame % len(self.atoms)
        self.redraw()

    frame = property(get_frame, set_frame, doc="Get or set the current frame")

    def first(self):
        """
        Show the first frame (frame 0).
        """

        self.frame = 0

    def last(self):
        """
        Show the last frame, i.e. len(self.atoms)-1
        """

        self.frame = len(self.atoms)-1

    def forward(self, delta=None):
        """
        Move forward by `delta` frames (default value is self.delta).
        """

        delta = delta or self.delta
        self.frame = ((self.frame + delta) % len(self.atoms))

    def backward(self, delta=None):
        """
        Move backward by `delta` frames (default values is self.delta).
        """
        delta = delta or self.delta
        self.frame = (self.frame - delta) % len(self.atoms)


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
        Move the crystal origin to fractional coordinates `s`

        For example, use ``s=[0.5, 0.5, 0.5]`` to shift by half the cell along
        the :math:`\mathbf{a}`, :math:`\mathbf{b}` and :math:`\mathbf{c}`
        lattice vectors.
        """
        sa, sb, sc = s
        self.run_command("xtal_origin_goto %f %f %f" % (sa, sb, sc))

    def find_atom(self, i):
        """
        Set the anchor to the atom with index `i`.
        """
        if get_fortran_indexing(): i = i-1
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

    def draw_arrows(self, property, scale_factor=0.0, head_height=0.1,
                    head_width=0.05, up=(0.0,1.0,0.0)):
        """
        Draw arrows on each atom, based on a vector property

        Parameters
        ----------

        property : string
           Name of the array to use for arrow vectors.
           Use ``None`` to turn off previous arrows.
        scale_factor : float
           Override length of arrows. 1 unit = 1 Angstrom; default
           value of 0.0 means autoscale.
        head_height : float
           Specify height of arrow heads in Angstrom. 
        head_width : float
        up : 3-vector (tuple, list or array)
           Specify the plane in which the arrow heads are
           drawn. Arrows are drawn in the plane which is common
           to their direction and this vector.
           Default is ``[0.,1.,0.]``.
        """
        
        if property is None:
            self.run_command('draw_arrows off')
        else:
            property = self._property_hook(self.gcat(), property)
            up1, up2, up3 = up
            self.redraw() # ensure property is available
            self.run_command('draw_arrows %s %f %f %f %f %f %f' %
                             (str(property), scale_factor, head_height, head_width, up1, up2, up3))

    def wait(self):
        """Sleep until this AtomEye viewer has finished processing all queued events."""
        if not self.is_alive: 
            raise RuntimeError('is_alive is False')
        self._atomeye.wait(self._window_id)


    def get_visible(self):
        """Return list of indices of atoms currently visible in this viewer."""
        indices = self._atomeye.get_visible()
        at = self.gcat()
        if np.any(indices > len(at)):
            # atoms duplicated due to narrow cell (n->small_cell_err_handler == 1)
            indices = list(set([idx % len(at) for idx in indices ]))
        if get_fortran_indexing():
            indices = [idx+1 for idx in indices]
        return indices


    def get_size_pixels(self, state=None):
        """
        Return (width, height) in pixels of this viewer
        """
        if state is None:
            state = self.get_state()
        for command in state['commands']:
            if command.startswith('resize'):
                break
        else:
            raise ValueError('cannot find "resize" entry in state["commands"]')
        resize, width, height = command.split()
        width = int(width)
        height = int(height)
        return (width, height)


    def get_size_angstrom(self, state=None):
        """
        Return (width, height) in Angstrom of currently projected view

        Assumes object lies in plane z=0
        """
        if state is None:
            state = self.get_state()
        W, H = self.get_size_pixels(state)

        # camera position
        cx, cy, cz = np.array([float(f) for f in state['variables']['AX_3D->x'].split()])

        # conversion factor from view angle to pixels
        k = float(state['variables']['AX_3D->k'])

        w = abs(W*cz/k)
        h = abs(H*cz/k)
        
        return (w, h)
     

    def display(self):
        """
        Display snapshot from AtomEye session in IPython notebook
        """
        fd,fname = tempfile.mkstemp(suffix='.png')
        os.close(fd)
        self.capture(fname)
        self.wait()
        from IPython.display import Image, display
        display(Image(filename=fname))
        os.unlink(fname)
        

_viewer = None

def gcv():
    return _viewer

def scv(viewer):
    _viewer = viewer

def view(atoms, **kwargs):
    """
    Convenience wrapper to create/reuse a default `AtomEyeViewer`
    """
    global _viewer
    if _viewer is None:
        _viewer = AtomEyeViewer(atoms, **kwargs)
    else:
        _viewer.show(atoms, **kwargs)
    return _viewer

def redraw():
    """
    Redraw current AtomEye window, keeping Atoms and settings the same.
    """
    gcv().redraw()

def run_command(command):
    """
    Run a command in current AtomEye thread.

    The command is queued for later execution, unless :attr:`block` is True.

    Parameters
    ----------

    command : string
       The command to pass to AtomEye
    """
    gcv().run_command(command)

def run_script(script):
    """
    Run commands from the file script, in a blocking fashion.
    """
    gcv().run_script(script)

def close():
    """
    Close the current viewer window.
    """
    gcv().close()

def setp(self, key, value):
    """
    Run the AtomEye command "set key value"
    """
    gcv().setp(key, value)

def save_script(filename):
    """
    Save AtomEye viewer settings to a file.
    """
    gcv().save(filename)

def toggle_coordination_coloring():
    """
    Turn on or off colouring by coordination number (key "k")
    """
    gcv().toggle_coordination_coloring()

def translate(axis, delta):
    """
    Translate system along `axis` by an amount `delta` (key "Ctrl+left/right/up/down")
    """
    gcv().translate(axis, delta)

def shift_xtal(axis, delta):
    """
    Shift crystal within periodic boundaries along `axis` by `delta` (key "Shift+left/right/up/down").
    """
    gcv().shift_xtal(axis, delta)

def rotate(axis, theta):
    """
    Rotate around `axis` by angle `theta`.
    """
    gcv().rotate(axis, theta)

def advance(delta):
    """
    Move the camera forward by `delta`.
    """
    gcv().advance(delta)

def shift_cutting_plane(delta):
    """
    Shift the current cutting plane by an amount `delta`.
    """
    gcv().shift_cutting_plane(delta)

def change_bgcolor(color):
    """
    Change the viewer background colour to `color`, which should be a RGB tuple with three floats in range 0..1.
    """
    gcv().change_bgcolor(color)

def change_atom_r_ratio(delta):
    """
    Change the size of the balls used to draw the atoms by `delta`.
    """
    gcv().change_atom_r_ratio(delta)

def change_bond_radius(delta):
    """
    Change the radius of the cylinders used the draw bonds by `delta`.
    """
    gcv().change_bond_radius(delta)

def change_view_angle_amplification(delta):
    """
    Change the amplification of the view angle by `delta`.
    """
    gcv().change_view_angle_amplification(delta)

def toggle_parallel_projection():
    """
    Toggle between parallel and perspective projections.
    """
    gcv().toggle_parallel_projection()

def toggle_bond_mode():
    """
    Turn on or off bonds.
    """
    gcv().toggle_bond_mode()

def toggle_small_cell_mode():
    """
    Toggle between two different behaviours for when cell is smaller than r_cut/2:
     1. clip cell - some neigbours may be lost (default)
     2. replicate cell along narrow directions
    """
    gcv().toggle_small_cell_mode()

def normal_coloring():
    """
    Return to normal colouring of the atoms (key "o").
    """
    gcv().normal_coloring()

def aux_property_coloring(auxprop):
    """
    Colour the currently viewed atoms according to `auxprop`.

    Overloaded to allow 
    See :ref:`qlab_atom_coloring` for more details and examples.

    Parameters
    ----------
    auxprop : str, array_like, int or list
       Values to use to colour the atoms. Should be either the
       name of a scalar field entry in :attr:`~.Atoms.properties`
       (or equivalently, :attr:`~Atoms.arrays`) such as ``"charge"``,
       a float, int or bool array of shape ``(len(gcat()),)``, or an
       atom index or list of atom indices to highlight particular atoms.
    """
    gcv().aux_property_coloring(auxprop)

def central_symmetry_coloring():
    """
    Colour atoms by centro-symmetry parameter.
    """
    gcv().central_symmetry_coloring()

def change_aux_property_threshold(lower, upper):
    """
    Change the lower and upper aux property thresholds.
    """
    gcv().change_aux_property_threshold(lower, upper)

def reset_aux_property_thresholds():
    """
    Reset aux property thresholds to automatic values.
    """
    gcv().reset_aux_property_thresholds()

def toggle_aux_property_thresholds_saturation():
    """
    Toggle between saturated colouring and invisibility for values outside aux prop thresholds.
    """
    gcv().toggle_aux_property_thresholds_saturation()

def toggle_aux_property_thresholds_rigid():
    """
    Toggle between floating and rigid aux property thresholds when moving between frames
    """
    gcv().toggle_aux_property_thresholds_rigid()

def rcut_patch(sym1, sym2, value, absolute=False):
    """
    Change the cutoff distance for `sym1`--`sym2` bonds by `delta`.

    e.g. to increase cutoff for Si-Si bonds by 0.5 A use::

         viewer.rcut_patch('Si', 'Si', 0.5)

    With `absolute` set to True, `value` is used to set the
    absolute cutoff distance for `sym1`--`sym2` bonds, e.g.::

         viewer.rcut_patch('Si', 'Si', 2.50, True)
    """
    gcv().rcut_patch(sym1, sym2, value, absolute)

def select_gear(gear):
    """
    Change the AtomEye gear to `gear`

    Equivalent to pressing the one of the numeric keys 0..9
    """
    gcv().select_gear(gear)

def cutting_plane(n, d, s):
    """
    Create a new cutting plane with index `n`, normal `d`, and fractional displacement `s`.
    """
    gcv().cutting_plane(n, d, s)

def shift_cutting_plane_to_anchor(n):
    """
    Move the cutting plane with index `n` to the anchor
    """
    gcv().shift_cutting_plane_to_anchor(n)

def delete_cutting_plane(n):
    """
    Delete the cutting plane with index `n`
    """
    gcv().delete_cutting_plane(n)

def flip_cutting_plane(n):
    """
    Flip the cutting plane with index `n`
    """
    gcv().flip_cutting_plane(n)

def capture(filename, resolution=None):
    """
    Render the current view to image `filename`

    Format is determined from file extension: .png, .jpeg, or .eps.
    """
    gcv().capture(filename, resolution)

def change_wireframe_mode():
    """
    Change the display mode for the unit cell box.

    Equivalent to pressing the `i` key.
    """
    gcv().change_wireframe_mode()

def change_cutting_plane_wireframe_mode():
    """
    Change the display mode for cutting planes
    """
    gcv().change_cutting_plane_wireframe_mode()

def get_frame():
    """
    Get index of frame currently being viewed
    """
    return gcv().frame

def set_frame(frame):
    """
    Set current frame index to `frame`
    """
    gcv().frame = frame

def get_delta():
    """
    Get frame increment rate
    """
    return gcv().delta

def set_delta(delta):
    """
    Set frame increment rate
    """
    gcv().delta = delta

def first():
    """
    Show the first frame (frame 0).
    """
    gcv().first()

def last():
    """
    Show the last frame, i.e. len(gcv())-1
    """
    gcv().last()

def forward(delta=None):
    """
    Move forward by `delta` frames (default value is gcv().delta).
    """
    gcv().forward(delta)

def backward(delta=None):
    """
    Move backward by `delta` frames (default values is gcv().delta).
    """
    gcv().backward(delta)

def load_atom_color(filename):
    """
    Load atom colours from a .clr file.
    """
    gcv().load_atom_color(filename)

def load_aux(filename):
    """
    Load aux property values from a .aux file.
    """
    gcv().load_aux(filename)

def look_at_the_anchor():
    """
    Equivalent to pressing the `a` key
    """
    gcv().look_at_the_anchor()

def observer_goto():
    """
    Prompt for fractional position and move the observer there

    Equivalent to pressing the `g` key.
    """
    gcv().observer_goto()

def xtal_origin_goto(s):
    """
    Move the crystal origin to fractional coordinates `s`

    For example, use ``s=[0.5, 0.5, 0.5]`` to shift by half the cell along
    the :math:`\mathbf{a}`, :math:`\mathbf{b}` and :math:`\mathbf{c}`
    lattice vectors.
    """
    gcv().xtal_origin_goto(s)

def find_atom(i):
    """
    Set the anchor to the atom with index `i`.
    """
    gcv().find_atom(i)

def resize(width, height):
    """
    Resize the current window to `width` x `height` pixels.
    """
    gcv().resize(width, height)

def change_aux_colormap(n):
    """
    Select the `n`\ -th auxiliary property colourmap. 
    """
    gcv().change_aux_colormap(n)

def draw_arrows(property, scale_factor=0.0, head_height=0.1,
                head_width=0.05, up=(0.0,1.0,0.0)):
    """
    Draw arrows on each atom, based on a vector property

    Parameters
    ----------
    property : string
       Name of the array to use for arrow vectors.
       Use ``None`` to turn off previous arrows.
    scale_factor : float
       Override length of arrows. 1 unit = 1 Angstrom; default
       value of 0.0 means autoscale.
    head_height : float
       Specify height of arrow heads in Angstrom. 
    head_width : float
    up : 3-vector (tuple, list or array)
       Specify the plane in which the arrow heads are
       drawn. Arrows are drawn in the plane which is common
       to their direction and this vector.
       Default is ``[0.,1.,0.]``.
    """
    gcv().draw_arrows(property, scale_factor, head_height,
                      head_width, up)

def wait():
    """Sleep until current AtomEye viewer has finished processing all queued events."""
    gcv().wait()

def display():
    """
    Display snapshot from current viewer in Jupyter notebook
    """
    gcv().display()


