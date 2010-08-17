// HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// HQ X
// HQ X   quippy: Python interface to QUIP atomistic simulation library
// HQ X
// HQ X   Copyright James Kermode 2010
// HQ X
// HQ X   These portions of the source code are released under the GNU General
// HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
// HQ X
// HQ X   If you would like to license the source code under different terms,
// HQ X   please contact James Kermode, james.kermode@gmail.com
// HQ X
// HQ X   When using this software, please cite the following reference:
// HQ X
// HQ X   http://www.jrkermode.co.uk/quippy
// HQ X
// HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include <Python.h>
#include <pthread.h>
#include <signal.h>

#include "numpy/arrayobject.h"

#include <libatoms.h>
#include <atomeyelib.h>

static char atomeye_doc[] = 
"This module interfaces to AtomEye.";

static PyObject *on_click_atom_pyfunc = NULL;
static PyObject *on_advance_pyfunc = NULL;
static PyObject *on_close_pyfunc = NULL;
static PyObject *on_new_pyfunc = NULL;

static Atomeyelib_atoms atomeye_atoms;

static void on_click_atom(int iw, int atom)
{
  PyObject *arglist;
  PyGILState_STATE state;

  state = PyGILState_Ensure();
  if (on_click_atom_pyfunc != NULL) {
    arglist = Py_BuildValue("(i,i)", iw, atom);       
    PyEval_CallObject(on_click_atom_pyfunc, arglist); 
    Py_DECREF(arglist);                               
  }
  PyGILState_Release(state);
}

static void on_close(int iw)
{
  PyObject *arglist;
  PyGILState_STATE state;

  state = PyGILState_Ensure();
  if (on_close_pyfunc != NULL) {
    arglist = Py_BuildValue("(i)", iw);
    PyEval_CallObject(on_close_pyfunc, arglist);
    Py_DECREF(arglist);
  }
  PyGILState_Release(state);
}

static void on_advance(int iw, char *instr)
{
  PyObject *arglist;
  PyGILState_STATE state;

  state = PyGILState_Ensure();
  if (on_advance_pyfunc != NULL) {
    arglist = Py_BuildValue("(i,s)", iw, instr);
    PyEval_CallObject(on_advance_pyfunc, arglist);
    Py_DECREF(arglist);                           
  }
  PyGILState_Release(state);
}

static void on_new(int iw)
{
  PyObject *arglist;
  PyGILState_STATE state;

  state = PyGILState_Ensure();
  if (on_new_pyfunc != NULL) {
    arglist = Py_BuildValue("(i)", iw);       
    PyEval_CallObject(on_new_pyfunc, arglist);
    Py_DECREF(arglist);                       
  }
  PyGILState_Release(state);
}

static int update_atoms_structure(PyObject *pyat)
{
  int i,j;
  PyObject *n = NULL, *lattice = NULL, *properties = NULL, *fpointer = NULL;

  /* atoms.n - int */
  if ((n = PyObject_GetAttrString(pyat, "n")) == NULL) {
    PyErr_SetString(PyExc_AttributeError, "Missing attribute atoms.n");
    goto fail;
  }
  if (!PyInt_Check(n)) {
    PyErr_SetString(PyExc_TypeError, "atoms.n must be an integer");
    goto fail;
  }
  atomeye_atoms.n_atom = (int)PyInt_AsLong(n);

  /* atoms.properties._fpointer */
  if ((properties = PyObject_GetAttrString(pyat, "properties")) == NULL) {
    PyErr_SetString(PyExc_AttributeError, "Missing attribute atoms.properties");
    goto fail;
  }
  if ((fpointer = PyObject_GetAttrString(properties, "_fpointer")) == NULL) {
    PyErr_SetString(PyExc_AttributeError, "Missing attribute atoms.properties._fpointer");
    goto fail;
  }
  if (PyArray_NDIM(fpointer) != 1) {
    PyErr_SetString(PyExc_ValueError, "atoms.properties._fpointer must have 1 dimension");
    goto fail;
  }
  if (PyArray_DIM(fpointer, 0) != SIZEOF_FORTRAN_T) {
    PyErr_SetString(PyExc_ValueError, "atoms.properties._fpointer must have shape (SIZEOF_FORTRAN_T,)");
    goto fail;
  }
  if (PyArray_TYPE(fpointer) != NPY_INT) {
    PyErr_SetString(PyExc_ValueError, "atoms.properties._fpointer must have type int");
    goto fail;
  }
  for (i=0; i < SIZEOF_FORTRAN_T; i++)
    atomeye_atoms.properties[i] = ((int*)PyArray_DATA(fpointer))[i];

  /* atoms.lattice */
  if ((lattice = PyObject_GetAttrString(pyat, "lattice")) == NULL) {
    PyErr_SetString(PyExc_AttributeError, "Missing atoms.lattice");
    goto fail;
  }
  if (PyArray_NDIM(lattice) != 2) {
    PyErr_SetString(PyExc_ValueError, "atoms.lattice must have 2 dimensions");
    goto fail;
  }
  if (PyArray_DIM(lattice,0) != 3 || PyArray_DIM(lattice,1) != 3) {
    PyErr_SetString(PyExc_ValueError, "atoms.lattice must have shape (3,3)");
    goto fail;
  }
  if (PyArray_TYPE(lattice) != NPY_DOUBLE) {
    PyErr_SetString(PyExc_ValueError, "atoms.lattice must have type double");
    goto fail;
  }
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      atomeye_atoms.lattice[i][j] = *(double *)PyArray_GETPTR2(lattice, i, j);

  Py_DECREF(n);                       
  Py_DECREF(properties);
  Py_DECREF(fpointer);
  Py_DECREF(lattice);
  return 1;

 fail:
  Py_DECREF(n);                       
  Py_DECREF(properties);
  Py_DECREF(fpointer);
  Py_DECREF(lattice);
  return 0;
}

static char atomeye_open_window_doc[]=
  "iw = _atomeye.open_window(copy=-1,atoms=None,nowindow=) -- open a new AtomEye window";

static PyObject*
atomeye_open_window(PyObject *self, PyObject *args)
{
  int icopy = -1, iw, argc;
  char outstr[255];
  char *argv[3];
  PyObject *pyat = NULL;
  static int atomeye_initialised = 0;
  int nowindow = 0;

  if (!PyArg_ParseTuple(args, "|iOi", &icopy, &pyat, &nowindow))
    return NULL;

  if (!atomeye_initialised) {
    argv[0] = (char *)malloc(20);
    argv[1] = (char *)malloc(20);
    strcpy(argv[0], "A");
    strcpy(argv[1], "-nostdin");
    argc = 2;
    if (nowindow) {
      strcpy(argv[2], "-nowindow");
      argc = 3;
    }
  
    if (pyat != NULL && pyat != Py_None) {
      if (!update_atoms_structure(pyat)) return NULL;
      atomeyelib_init(argc, argv, &atomeye_atoms);
    } else
      atomeyelib_init(argc, argv, NULL);

    atomeyelib_set_handlers(&on_click_atom, &on_close, &on_advance, &on_new);

    free(argv[0]);
    free(argv[1]);

    atomeye_initialised = 1;
  }

  iw = atomeyelib_open_window(icopy);

  if (iw == -1) {
    sprintf(outstr, "Bad copy window id %d", icopy);
    PyErr_SetString(PyExc_RuntimeError, outstr);
    return NULL;
  }

  return PyInt_FromLong((long)iw);
}

static char atomeye_set_handlers_doc[]=
  "_atomeye.set_handlers(on_click, on_close, on_advance)";

static PyObject*
atomeye_set_handlers(PyObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, "OOOO", &on_click_atom_pyfunc, &on_close_pyfunc, &on_advance_pyfunc, &on_new_pyfunc))
    return NULL;

  if (!PyCallable_Check(on_click_atom_pyfunc)) {
    PyErr_SetString(PyExc_TypeError, "Need a callable object!");
    return NULL;
  }

  if (!PyCallable_Check(on_close_pyfunc)) {
    PyErr_SetString(PyExc_TypeError, "Need a callable object!");
    return NULL;
  }

  if (!PyCallable_Check(on_advance_pyfunc)) {
    PyErr_SetString(PyExc_TypeError, "Need a callable object!");
    return NULL;
  }

  if (!PyCallable_Check(on_new_pyfunc)) {
    PyErr_SetString(PyExc_TypeError, "Need a callable object!");
    return NULL;
  }

  Py_INCREF(on_click_atom_pyfunc);
  Py_INCREF(on_advance_pyfunc);
  Py_INCREF(on_close_pyfunc);
  Py_INCREF(on_new_pyfunc);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static char atomeye_redraw_doc[]=
  "_atomeye.redraw(iw) -- redraw window";

static PyObject*
atomeye_redraw(PyObject *self, PyObject *args)
{
  int iw;
  char outstr[255];

  if (!PyArg_ParseTuple(args, "i", &iw))
    return NULL;

  if (!atomeyelib_queueevent(iw, ATOMEYELIB_REDRAW, "", NULL, outstr)) {
    PyErr_SetString(PyExc_RuntimeError, outstr);
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static char atomeye_run_command_doc[]=
  "_atomeye.run_command(iw, instr) -- run an AtomEye command";

static PyObject*
atomeye_run_command(PyObject *self, PyObject *args)
{
  int iw;
  char *instr;
  char outstr[255];
  
  if (!PyArg_ParseTuple(args, "is", &iw, &instr))
    return NULL;

  if (!atomeyelib_queueevent(iw, ATOMEYELIB_RUN_COMMAND, instr, NULL, outstr)) {
    PyErr_SetString(PyExc_RuntimeError, outstr);
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static char atomeye_load_atoms_doc[]=
  "_atomeye.load_atoms(iw, title, atoms) -- load atoms into AtomEye";

static PyObject*
atomeye_load_atoms(PyObject *self, PyObject *args)
{
  int iw;
  char *title;
  PyObject *pyat;
  char outstr[255];
  
  if (!PyArg_ParseTuple(args, "isO", &iw, &title, &pyat))
    return NULL;
  
  if (!update_atoms_structure(pyat)) return NULL;

  if (!atomeyelib_queueevent(iw, ATOMEYELIB_LOAD_ATOMS, title, &atomeye_atoms, outstr)) {
    PyErr_SetString(PyExc_RuntimeError, outstr);    
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static char atomeye_set_title_doc[] =
  "_atomeye.set_title(iw, title) -- set window title";

static PyObject*
atomeye_set_title(PyObject *self, PyObject *args)
{
  int iw;
  char *title;

  if (!PyArg_ParseTuple(args, "is", &iw, &title))
    return NULL;

  atomeyelib_set_title(iw, title);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static char atomeye_wait_doc[] =
  "_atomeye.wait(iw) -- wait for window `iw` to complete all queued events";

static PyObject*
atomeye_wait(PyObject *self, PyObject *args)
{
  int iw;

  if (!PyArg_ParseTuple(args, "i", &iw))
    return NULL;
  
  atomeyelib_wait(iw);
  return Py_None;
}


static PyMethodDef atomeye_methods[] = {
  {"open_window", atomeye_open_window, METH_VARARGS, atomeye_open_window_doc},
  {"set_handlers", atomeye_set_handlers, METH_VARARGS, atomeye_set_handlers_doc},
  {"redraw", atomeye_redraw, METH_VARARGS, atomeye_redraw_doc},
  {"run_command", atomeye_run_command, METH_VARARGS, atomeye_run_command_doc},
  {"load_atoms", atomeye_load_atoms, METH_VARARGS, atomeye_load_atoms_doc},
  {"set_title", atomeye_set_title, METH_VARARGS, atomeye_set_title_doc},
  {"wait", atomeye_wait, METH_VARARGS, atomeye_wait_doc},
  {NULL, NULL}
};

PyMODINIT_FUNC
init_atomeye(void)
{
  Py_InitModule3("_atomeye", atomeye_methods, atomeye_doc);
  import_array();
}
