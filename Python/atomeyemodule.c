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
static PyObject *on_redraw_pyfunc = NULL;

static int update_atoms_structure(int n_atom, PyObject *cell, PyObject *arrays, Atomeyelib_atoms *atoms);

static void on_click_atom(int mod_id, int iw, int atom)
{
  PyObject *arglist;
  PyGILState_STATE state;

  state = PyGILState_Ensure();
  if (on_click_atom_pyfunc != NULL) {
    arglist = Py_BuildValue("(i,i,i)", mod_id, iw, atom);       
    PyEval_CallObject(on_click_atom_pyfunc, arglist); 
    Py_DECREF(arglist);                               
  }
  PyGILState_Release(state);
}

static void on_close(int mod_id, int iw)
{
  PyObject *arglist;
  PyGILState_STATE state;

  state = PyGILState_Ensure();
  if (on_close_pyfunc != NULL) {
    arglist = Py_BuildValue("(i,i)", mod_id, iw);
    PyEval_CallObject(on_close_pyfunc, arglist);
    Py_DECREF(arglist);
  }
  PyGILState_Release(state);
}

static void on_advance(int mod_id, int iw, char *instr)
{
  PyObject *arglist;
  PyGILState_STATE state;

  state = PyGILState_Ensure();
  if (on_advance_pyfunc != NULL) {
    arglist = Py_BuildValue("(i,i,s)", mod_id, iw, instr);
    PyEval_CallObject(on_advance_pyfunc, arglist);
    Py_DECREF(arglist);                           
  }
  PyGILState_Release(state);
}

static void on_new(int mod_id, int iw)
{
  PyObject *arglist;
  PyGILState_STATE state;

  state = PyGILState_Ensure();
  if (on_new_pyfunc != NULL) {
    arglist = Py_BuildValue("(i,i)", mod_id, iw);
    PyEval_CallObject(on_new_pyfunc, arglist);
    Py_DECREF(arglist);                       
  }
  PyGILState_Release(state);
}

static int on_redraw(int mod_id, int iw, Atomeyelib_atoms *atoms)
{
  PyObject *arglist, *result;
  PyGILState_STATE state;
  int redraw, nat;
  char *mytitle;
  PyObject *cell, *arrays;

  redraw = 0;
  state = PyGILState_Ensure();
  if (on_redraw_pyfunc != NULL) {
    arglist = Py_BuildValue("(i,i)", mod_id, iw);
    result = PyEval_CallObject(on_redraw_pyfunc, arglist);
    Py_DECREF(arglist);                       

    if (result == NULL) {
      fprintf(stderr, "WARNING: on_redraw(mod_id=%d, iw=%d) returned NULL\n", mod_id, iw);
      return 0;
    }

    if (!PyArg_ParseTuple(result, "isiOO", &redraw, &mytitle, &nat, &cell, &arrays)) {
      PyErr_PrintEx(0);
      return 0;
    }
    strcpy(atoms->title, mytitle);
    if (!redraw) {
      // nothing to be done
      Py_DECREF(result);
      PyGILState_Release(state);
      return 0;
    }

    if (!update_atoms_structure(nat, cell, arrays, atoms)) {
      PyErr_PrintEx(0);
      Py_DECREF(result);
      PyGILState_Release(state);
      return 0;
    }
    Py_DECREF(result);
  }
  PyGILState_Release(state);
  return redraw;
}

static int update_atoms_structure(int n_atom, PyObject *cell, PyObject *arrays, Atomeyelib_atoms *atoms)
{
  int i,j, type, shape[2], error = 0;
  PyObject *fpointer = NULL, *pykey = NULL, *array = NULL, *iterator = NULL;
  char *key;
  void *data;

  atoms->n_atom = n_atom;

  /* cell - 3x3 array of lattice vector as rows */
  if (PyArray_NDIM(cell) != 2) {
    PyErr_SetString(PyExc_ValueError, "cell must have 2 dimensions");
    goto fail;
  }
  if (PyArray_DIM(cell,0) != 3 || PyArray_DIM(cell,1) != 3) {
    PyErr_SetString(PyExc_ValueError, "cell must have shape (3,3)");
    goto fail;
  }
  if (PyArray_TYPE(cell) != NPY_DOUBLE) {
    PyErr_SetString(PyExc_ValueError, "cell must have type double");
    goto fail;
  }
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      atoms->lattice[i][j] = *(double *)PyArray_GETPTR2(cell, i, j);

  /* Test for arrays._fpointer */
  if (!PyObject_HasAttrString(arrays, "_fpointer")) {
    if (!PyMapping_Check(arrays)) {
      PyErr_SetString(PyExc_ValueError, "arrays must support mapping protocol");
      goto fail;
    }

    /* construct new Fortran dictionary, copying keys and values from arrays */
    dictionary_initialise(atoms->properties);
    atoms->allocated = 1;
    
    /* for key in arrays */
    iterator = PyObject_GetIter(arrays);
    while ((pykey = PyIter_Next(iterator))) {
      if (!PyString_Check(pykey)) {
	PyErr_SetString(PyExc_ValueError, "All keys in arrays must be strings");
	goto fail;
      }
      key = PyString_AsString(pykey);
      /* array = arrays[key] */
      if ((array = PyMapping_GetItemString(arrays, key)) == NULL) { 
	PyErr_Format(PyExc_ValueError, "Error accessing key %s in arrays", key);
	goto fail;
      }
      if (!(PyArray_TYPE(array) == NPY_FLOAT32 || PyArray_TYPE(array) == NPY_FLOAT64 ||
	    PyArray_TYPE(array) == NPY_INT32   || PyArray_TYPE(array) == NPY_INT64)) {
	Py_DECREF(pykey);
	Py_DECREF(array);
	continue;
      }
      
      if (PyArray_NDIM(array) == 1) {
	if (PyArray_DIM(array, 0) != n_atom) {
	  PyErr_Format(PyExc_ValueError, "Array \"%s\" in arrays must have size n_atoms (%d)", key, n_atom);
	  goto fail;
	}
	
	if (PyArray_TYPE(array) == NPY_FLOAT32 || PyArray_TYPE(array) == NPY_FLOAT64) 
	  type = T_REAL_A;
	else
	  type = T_INTEGER_A;
	shape[0] = n_atom;
	dictionary_add_key(atoms->properties, key, &type, shape, &data, &error, strlen(key));
	if (error != 0) {
	  PyErr_Format(PyExc_ValueError, "Error adding key \"%s\" shape (%d,%d)", key, shape[0], shape[1]);
	  goto fail;
	}
	if (type == T_REAL_A) {
	  for (i=0; i<n_atom; i++)
	    REAL_A(data, i) = *(double *)PyArray_GETPTR1(array, i);
	} else {
	  for (i=0; i<n_atom; i++)
	    INTEGER_A(data, i) = *(int *)PyArray_GETPTR1(array, i);
	}
      } else if (PyArray_NDIM(array) == 2) {
	if (PyArray_DIM(array, 0) != n_atom) {
	  PyErr_Format(PyExc_ValueError, "Array \"%s\" in arrays must have shape (n_atom,n_cols)", key);
	  goto fail;
	}	

	if (PyArray_TYPE(array) == NPY_FLOAT32 || PyArray_TYPE(array) == NPY_FLOAT64) 
	  type = T_REAL_A2;
	else
	  type = T_INTEGER_A2;
	shape[0] = PyArray_DIM(array, 1);
	shape[1] = n_atom;
	dictionary_add_key(atoms->properties, key, &type, shape, &data, &error, strlen(key));	
	if (error != 0) {
	  PyErr_Format(PyExc_ValueError, "Error adding key \"%s\" shape (%d,%d)", key, shape[0], shape[1]);
	  goto fail;
	}
	if (type == T_REAL_A2) {
	  for (j=0; j<shape[1]; j++)
	    for (i=0; i<shape[0]; i++)
	      REAL_A2(data, shape, i, j) = *(double *)PyArray_GETPTR2(array, j, i);
	} else {
	    for (j=0; j<shape[1]; j++)
	      for (i=0; i<shape[0]; i++)
		INTEGER_A2(data, shape, i, j) = *(int *)PyArray_GETPTR2(array, j, i);
	}
      } else {
	PyErr_SetString(PyExc_ValueError, "Array \"%s\" in arrays must have either 1 or 2 dimensions");
	goto fail;
      }

      Py_DECREF(pykey);
      Py_DECREF(array);
    }
    
    Py_DECREF(iterator);

  } else {
    /* arrays has _fpointer attribute, so it's already a Fortran dictionary of arrays  */
    fpointer = PyObject_GetAttrString(arrays, "_fpointer");
    atoms->allocated = 0; // we didn't allocate dictionary, so shouldn't deallocate it
    if (PyArray_NDIM(fpointer) != 1) {
      PyErr_SetString(PyExc_ValueError, "arrays._fpointer must have 1 dimension");
      goto fail;
    }
    if (PyArray_DIM(fpointer, 0) != SIZEOF_FORTRAN_T) {
      PyErr_SetString(PyExc_ValueError, "arrays._fpointer must have shape (SIZEOF_FORTRAN_T,)");
      goto fail;
    }
    if (PyArray_TYPE(fpointer) != NPY_INT) {
      PyErr_SetString(PyExc_ValueError, "arrays._fpointer must have type int");
      goto fail;
    }
    for (i=0; i < SIZEOF_FORTRAN_T; i++)
      atoms->properties[i] = ((int*)PyArray_DATA(fpointer))[i];

    Py_DECREF(fpointer);
  }

  return 1;

 fail:
  Py_XDECREF(fpointer);
  Py_XDECREF(pykey);
  Py_XDECREF(array);
  Py_XDECREF(iterator);
  return 0;
}

static char atomeye_open_window_doc[]=
  "iw = _atomeye.open_window(mod_id=0, copy=-1, nowindow=0) -- open a new AtomEye window";

static PyObject*
atomeye_open_window(PyObject *self, PyObject *args)
{
  int mod_id = 0, icopy = -1, iw, argc;
  char outstr[255];
  char *argv[3];
  Atomeyelib_atoms atoms;
  static int atomeye_initialised = 0;
  int nowindow = 0;

  if (!PyArg_ParseTuple(args, "|iii", &mod_id, &icopy, &nowindow))
    return NULL;

  if (!atomeye_initialised) {
    atoms.allocated = 0;

    argv[0] = (char *)malloc(20);
    argv[1] = (char *)malloc(20);
    argv[2] = (char *)malloc(20);
    strcpy(argv[0], "A");
    strcpy(argv[1], "-nostdin");
    argc = 2;
    if (nowindow) {
      strcpy(argv[2], "-nowindow");
      argc = 3;
    }
  
    if (on_redraw(mod_id, 0, &atoms)) {
      atomeyelib_init(argc, argv, &atoms);
    } else
      atomeyelib_init(argc, argv, NULL);

    atomeyelib_set_handlers(&on_click_atom, &on_close, &on_advance, &on_new, &on_redraw);

    if (atoms.allocated)
      dictionary_finalise(atoms.properties);
    free(argv[0]);
    free(argv[1]);
    free(argv[2]);

    atomeye_initialised = 1;
  }

  iw = atomeyelib_open_window(mod_id, icopy);

  if (iw == -1) {
    sprintf(outstr, "Bad copy window id %d", icopy);
    PyErr_SetString(PyExc_RuntimeError, outstr);
    return NULL;
  } else if (iw == -2) {
    sprintf(outstr, "Cannot allocate sufficient stack space for AtomEye thread");
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
  if (!PyArg_ParseTuple(args, "OOOOO", &on_click_atom_pyfunc, &on_close_pyfunc, 
			&on_advance_pyfunc, &on_new_pyfunc, &on_redraw_pyfunc))
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

  if (!PyCallable_Check(on_redraw_pyfunc)) {
    PyErr_SetString(PyExc_TypeError, "Need a callable object!");
    return NULL;
  }

  Py_INCREF(on_click_atom_pyfunc);
  Py_INCREF(on_advance_pyfunc);
  Py_INCREF(on_close_pyfunc);
  Py_INCREF(on_new_pyfunc);
  Py_INCREF(on_redraw_pyfunc);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static char atomeye_redraw_doc[]=
  "_atomeye.redraw(iw) -- update atoms and redraw window";

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
  
  Py_BEGIN_ALLOW_THREADS
  atomeyelib_wait(iw);
  Py_END_ALLOW_THREADS
  return Py_None;
}

static char atomeye_get_visible_doc[] =
  "_atomeye.get_visible() -- return list of all visible atoms";

static PyObject*
atomeye_get_visible(PyObject *self, PyObject *args)
{
  int n_shown, *idx;
  PyObject *list = NULL, *entry = NULL;
  int i;
  
  atomeyelib_get_visible(&n_shown, &idx); //allocates pointer *idx
  list = PyList_New(0);
  for (i = 0; i < n_shown; i++) {
    entry = PyInt_FromLong(idx[i]);
    PyList_Append(list, entry);
    Py_DECREF(entry);
  }
  free(idx); // free temporary memory
  return list;
}

static PyMethodDef atomeye_methods[] = {
  {"open_window", atomeye_open_window, METH_VARARGS, atomeye_open_window_doc},
  {"set_handlers", atomeye_set_handlers, METH_VARARGS, atomeye_set_handlers_doc},
  {"redraw", atomeye_redraw, METH_VARARGS, atomeye_redraw_doc},
  {"run_command", atomeye_run_command, METH_VARARGS, atomeye_run_command_doc},
  {"set_title", atomeye_set_title, METH_VARARGS, atomeye_set_title_doc},
  {"wait", atomeye_wait, METH_VARARGS, atomeye_wait_doc},
  {"get_visible", atomeye_get_visible, METH_VARARGS, atomeye_get_visible_doc},
  {NULL, NULL}
};

PyMODINIT_FUNC
init_atomeye(void)
{
  Py_InitModule3("_atomeye", atomeye_methods, atomeye_doc);
  import_array();
}
