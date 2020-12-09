/* Copyright (C) 2020 Atsushi Togo */
/* All rights reserved. */

/* This file is part of kspclib. */

/* Redistribution and use in source and binary forms, with or without */
/* modification, are permitted provided that the following conditions */
/* are met: */

/* * Redistributions of source code must retain the above copyright */
/*   notice, this list of conditions and the following disclaimer. */

/* * Redistributions in binary form must reproduce the above copyright */
/*   notice, this list of conditions and the following disclaimer in */
/*   the documentation and/or other materials provided with the */
/*   distribution. */

/* * Neither the name of the kspclib project nor the names of its */
/*   contributors may be used to endorse or promote products derived */
/*   from this software without specific prior written permission. */

/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS */
/* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT */
/* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS */
/* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE */
/* COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, */
/* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, */
/* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; */
/* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER */
/* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT */
/* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN */
/* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE */
/* POSSIBILITY OF SUCH DAMAGE. */

#include <Python.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include <kspclib.h>

#if (PY_MAJOR_VERSION < 3) && (PY_MINOR_VERSION < 6)
#define PYUNICODE_FROMSTRING PyString_FromString
#else
#define PYUNICODE_FROMSTRING PyUnicode_FromString
#endif

static PyObject * py_version(PyObject *self, PyObject *args);
static PyObject * py_all_grid_addresses(PyObject *self, PyObject *args);
static PyObject * py_grid_point_double_mesh(PyObject *self, PyObject *args);
static PyObject * py_grid_address_double_mesh(PyObject *self, PyObject *args);
static PyObject * py_thm_relative_grid_addresses(PyObject *self, PyObject *args);
static PyObject * py_thm_integration_weight(PyObject *self, PyObject *args);
static PyObject * py_snf3x3(PyObject *self, PyObject *args);
static PyObject * py_snf_transform_rotations(PyObject *self, PyObject *args);

struct module_state {
  PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

static PyObject *
error_out(PyObject *m) {
  struct module_state *st = GETSTATE(m);
  PyErr_SetString(st->error, "something bad happened");
  return NULL;
}

static PyMethodDef _kspclib_methods[] = {
  {"error_out", (PyCFunction)error_out, METH_NOARGS, NULL},
  {"version", py_version, METH_VARARGS, "kspclib version"},
  {"all_grid_addresses", py_all_grid_addresses, METH_VARARGS,
   "Return all grid addresses for mesh"},
  {"grid_point_double_mesh", py_grid_point_double_mesh, METH_VARARGS,
   "Return grid point index of grid address of mesh"},
  {"grid_address_double_mesh", py_grid_address_double_mesh, METH_VARARGS,
   "Convert grid address plus shift to double-grid address"},
  {"thm_relative_grid_addresses", py_thm_relative_grid_addresses, METH_VARARGS,
   "Return relative grid addresses of 24 tetrahedra"},
  {"thm_integration_weight", py_thm_integration_weight, METH_VARARGS,
   "Return integration weight of tetrahedron method for a grid point"},
  {"snf3x3", py_snf3x3, METH_VARARGS,
   "Return D, P, Q of Smith normal form"},
  {"snf_transform_rotations", py_snf_transform_rotations, METH_VARARGS,
   "Transform rotations by SNF of grid generation matrix"},
  {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3

static int _kspclib_traverse(PyObject *m, visitproc visit, void *arg) {
  Py_VISIT(GETSTATE(m)->error);
  return 0;
}

static int _kspclib_clear(PyObject *m) {
  Py_CLEAR(GETSTATE(m)->error);
  return 0;
}

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_kspclib",
  NULL,
  sizeof(struct module_state),
  _kspclib_methods,
  NULL,
  _kspclib_traverse,
  _kspclib_clear,
  NULL
};

#define INITERROR return NULL

PyObject *
PyInit__kspclib(void)

#else
#define INITERROR return

  void
  init_kspclib(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
  PyObject *module = PyModule_Create(&moduledef);
#else
  PyObject *module = Py_InitModule("_kspclib", _kspclib_methods);
#endif

  if (module == NULL)
    INITERROR;
  struct module_state *st = GETSTATE(module);

  st->error = PyErr_NewException("_kspclib.Error", NULL, NULL);
  if (st->error == NULL) {
    Py_DECREF(module);
    INITERROR;
  }

#if PY_MAJOR_VERSION >= 3
  return module;
#endif
}

static PyObject * py_version(PyObject *self, PyObject *args)
{
  PyObject *array;
  int i;
  int version[3];

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  version[0] = ksp_get_major_version();
  version[1] = ksp_get_minor_version();
  version[2] = ksp_get_micro_version();

  array = PyList_New(3);
  for (i = 0; i < 3; i++) {
    PyList_SetItem(array, i, PyLong_FromLong((long)version[i]));
  }

  return array;
}

static PyObject * py_all_grid_addresses(PyObject *self, PyObject *args)
{
  PyArrayObject* py_grid_address;
  PyArrayObject* py_mesh;

  int (*grid_address)[3];
  int* mesh;

  if (!PyArg_ParseTuple(args, "OO",
                        &py_grid_address,
                        &py_mesh)) {
    return NULL;
  }

  grid_address = (int(*)[3])PyArray_DATA(py_grid_address);
  mesh = (int*)PyArray_DATA(py_mesh);

  ksp_get_all_grid_addresses(grid_address, mesh);

  Py_RETURN_NONE;
}

static PyObject * py_grid_point_double_mesh(PyObject *self, PyObject *args)
{
  PyArrayObject* py_address_double;
  PyArrayObject* py_mesh;

  int* address_double;
  int* mesh;
  size_t grid_point;


  if (!PyArg_ParseTuple(args, "OO",
                        &py_address_double,
                        &py_mesh)) {
    return NULL;
  }

  address_double = (int*)PyArray_DATA(py_address_double);
  mesh = (int*)PyArray_DATA(py_mesh);

  grid_point = ksp_get_grid_point_double_mesh(address_double, mesh);

  return PyLong_FromSize_t(grid_point);

}

static PyObject * py_grid_address_double_mesh(PyObject *self, PyObject *args)
{
  PyArrayObject* py_address_double;
  PyArrayObject* py_address;
  PyArrayObject* py_is_shift;
  PyArrayObject* py_mesh;

  int* address_double;
  int* address;
  int* is_shift;
  int* mesh;

  if (!PyArg_ParseTuple(args, "OOOO",
                        &py_address_double,
                        &py_address,
                        &py_mesh,
                        &py_is_shift)) {
    return NULL;
  }

  address_double = (int*)PyArray_DATA(py_address_double);
  address = (int*)PyArray_DATA(py_address);
  is_shift = (int*)PyArray_DATA(py_is_shift);
  mesh = (int*)PyArray_DATA(py_mesh);

  ksp_get_grid_address_double_mesh(address_double,
                                   address,
                                   mesh,
                                   is_shift);

  Py_RETURN_NONE;
}

static PyObject * py_thm_relative_grid_addresses(PyObject *self, PyObject *args)
{
  PyArrayObject* py_relative_grid_addresses;
  PyArrayObject* py_rec_lattice; /* column vectors */

  int (*relative_grid_addresses)[4][3];
  double (*rec_lattice)[3];

  if (!PyArg_ParseTuple(args, "OO",
                        &py_relative_grid_addresses,
                        &py_rec_lattice)) {
    return NULL;
  }

  relative_grid_addresses = (int(*)[4][3])PyArray_DATA(py_relative_grid_addresses);
  rec_lattice = (double(*)[3])PyArray_DATA(py_rec_lattice);

  ksp_get_thm_relative_grid_addresses(relative_grid_addresses,
                                      rec_lattice);

  Py_RETURN_NONE;
}

static PyObject * py_thm_integration_weight(PyObject *self, PyObject *args)
{
  double omega;
  PyArrayObject* py_tetrahedra_omegas;
  char* function;  /* I: delta function, J: Heaviside function */

  double (*tetrahedra_omegas)[4];  /* [24][4] */
  double iw;

  if (!PyArg_ParseTuple(args, "dOs",
                        &omega,
                        &py_tetrahedra_omegas,
                        &function)) {
    return NULL;
  }

  tetrahedra_omegas = (double(*)[4])PyArray_DATA(py_tetrahedra_omegas);

  iw = ksp_get_thm_integration_weight(omega,
                                      tetrahedra_omegas,
                                      function[0]);

  return PyFloat_FromDouble(iw);
}

static PyObject * py_snf3x3(PyObject *self, PyObject *args)
{
  PyArrayObject* py_DPQ;
  PyArrayObject* py_A;

  long (*DPQ)[3][3];  /* [3][3][3], left-most index gives D, P, Q. */
  long (*A)[3];  /* [3][3] */
  int succeeded;

  if (!PyArg_ParseTuple(args, "OO",
                        &py_DPQ,
                        &py_A)) {
    return NULL;
  }

  DPQ = (long(*)[3][3])PyArray_DATA(py_DPQ);
  A = (long(*)[3])PyArray_DATA(py_A);

  succeeded = ksp_get_snf3x3(DPQ[0], DPQ[1], DPQ[2], A);

  return PyBool_FromLong((long) succeeded);
}

static PyObject * py_snf_transform_rotations(PyObject *self, PyObject *args)
{
  PyArrayObject* py_transformed_rots;
  PyArrayObject* py_D;
  PyArrayObject* py_Q;
  PyArrayObject* py_rotations;

  long (*transformed_rots)[3][3];
  long (*D)[3];  /* [3][3] */
  long (*Q)[3];  /* [3][3] */
  int (*rotations)[3][3];
  int succeeded, num_rot;

  if (!PyArg_ParseTuple(args, "OOOO",
                        &py_transformed_rots,
                        &py_rotations,
                        &py_D,
                        &py_Q)) {
    return NULL;
  }

  transformed_rots = (long(*)[3][3])PyArray_DATA(py_transformed_rots);
  D = (long(*)[3])PyArray_DATA(py_D);
  Q = (long(*)[3])PyArray_DATA(py_Q);
  rotations = (int(*)[3][3])PyArray_DATA(py_rotations);
  num_rot = PyArray_DIMS(py_rotations)[0];

  succeeded = ksp_snf_transform_rotations(transformed_rots,
                                          rotations, num_rot, D, Q);

  return PyBool_FromLong((long) succeeded);
}
