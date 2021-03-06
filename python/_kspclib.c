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
static PyObject * py_double_grid_address(PyObject *self, PyObject *args);
static PyObject * py_double_grid_index(PyObject *self, PyObject *args);
static PyObject * py_thm_relative_grid_addresses(PyObject *self, PyObject *args);
static PyObject * py_thm_integration_weight(PyObject *self, PyObject *args);
static PyObject * py_snf3x3(PyObject *self, PyObject *args);
static PyObject * py_snf_transform_rotations(PyObject *self, PyObject *args);
static PyObject * py_all_grgrid_addresses(PyObject *self, PyObject *args);
static PyObject * py_double_grgrid_address(PyObject *self, PyObject *args);
static PyObject * py_grgrid_index(PyObject *self, PyObject *args);
static PyObject * py_double_grgrid_index(PyObject *self, PyObject *args);
static PyObject * py_grgrid_address_from_index(PyObject *self, PyObject *args);
static PyObject * py_rotate_grgrid_index(PyObject *self, PyObject *args);
static PyObject * py_ir_grgrid_map(PyObject *self, PyObject *args);
/* static PyObject * py_niggli_reduce(PyObject *self, PyObject *args); */
static PyObject * py_reciprocal_point_group(PyObject *self, PyObject *args);

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
   "Return all single-grid addresses"},
  {"double_grid_address", py_double_grid_address, METH_VARARGS,
   "Convert grid address plus shift to double-grid address"},
  {"double_grid_index", py_double_grid_index, METH_VARARGS,
   "Return grid point index of a double-grid address"},
  {"thm_relative_grid_addresses", py_thm_relative_grid_addresses, METH_VARARGS,
   "Return relative grid addresses of 24 tetrahedra"},
  {"thm_integration_weight", py_thm_integration_weight, METH_VARARGS,
   "Return integration weight of tetrahedron method for a grid point"},
  {"snf3x3", py_snf3x3, METH_VARARGS,
   "Return D, P, Q of Smith normal form"},
  {"snf_transform_rotations", py_snf_transform_rotations, METH_VARARGS,
   "Transform rotations by SNF of grid generation matrix"},
  {"all_grgrid_addresses", py_all_grgrid_addresses, METH_VARARGS,
   "Return all generalized regular grid addresses for diagonal elements of D"},
  {"double_grgrid_address", py_double_grgrid_address, METH_VARARGS,
   "Convert generalized-regular-grid address plus shift to double-grid address"},
  {"grgrid_index", py_grgrid_index, METH_VARARGS,
   "Return generalized-regular-grid point index of a single-grid address"},
  {"double_grgrid_index", py_double_grgrid_index, METH_VARARGS,
   "Return generalized-regular-grid point index of a double-grid address"},
  {"grgrid_address_from_index", py_grgrid_address_from_index, METH_VARARGS,
   "Return generalized-regular-grid address from its index"},
  {"rotate_grgrid_index", py_rotate_grgrid_index, METH_VARARGS,
   "Rotate a generalized-regular-grid point by index"},
  {"ir_grgrid_map", py_ir_grgrid_map, METH_VARARGS,
   "Find irreducible generalized-regular-grid points"},
  /* {"niggli_reduce", py_niggli_reduce, METH_VARARGS,
   *  "Perform Niggli reduction"}, */
  {"reciprocal_point_group", py_reciprocal_point_group, METH_VARARGS,
   "Return reciprocal point group"},
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
  long i;
  long version[3];

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

  long (*grid_address)[3];
  long* mesh;

  if (!PyArg_ParseTuple(args, "OO",
                        &py_grid_address,
                        &py_mesh)) {
    return NULL;
  }

  grid_address = (long(*)[3])PyArray_DATA(py_grid_address);
  mesh = (long*)PyArray_DATA(py_mesh);

  ksp_get_all_grid_addresses(grid_address, mesh);

  Py_RETURN_NONE;
}


static PyObject * py_double_grid_address(PyObject *self, PyObject *args)
{
  PyArrayObject* py_address_double;
  PyArrayObject* py_address;
  PyArrayObject* py_is_shift;
  PyArrayObject* py_mesh;

  long* address_double;
  long* address;
  long* is_shift;
  long* mesh;

  if (!PyArg_ParseTuple(args, "OOOO",
                        &py_address_double,
                        &py_address,
                        &py_mesh,
                        &py_is_shift)) {
    return NULL;
  }

  address_double = (long*)PyArray_DATA(py_address_double);
  address = (long*)PyArray_DATA(py_address);
  is_shift = (long*)PyArray_DATA(py_is_shift);
  mesh = (long*)PyArray_DATA(py_mesh);

  ksp_get_double_grid_address(address_double,
                              address,
                              mesh,
                              is_shift);

  Py_RETURN_NONE;
}


static PyObject * py_double_grid_index(PyObject *self, PyObject *args)
{
  PyArrayObject* py_address_double;
  PyArrayObject* py_mesh;

  long* address_double;
  long* mesh;
  long grid_index;


  if (!PyArg_ParseTuple(args, "OO",
                        &py_address_double,
                        &py_mesh)) {
    return NULL;
  }

  address_double = (long*)PyArray_DATA(py_address_double);
  mesh = (long*)PyArray_DATA(py_mesh);

  grid_index = ksp_get_double_grid_index(address_double, mesh);

  return PyLong_FromLong(grid_index);

}


static PyObject * py_thm_relative_grid_addresses(PyObject *self, PyObject *args)
{
  PyArrayObject* py_relative_grid_addresses;
  PyArrayObject* py_rec_lattice; /* column vectors */

  long (*relative_grid_addresses)[4][3];
  double (*rec_lattice)[3];

  if (!PyArg_ParseTuple(args, "OO",
                        &py_relative_grid_addresses,
                        &py_rec_lattice)) {
    return NULL;
  }

  relative_grid_addresses = (long(*)[4][3])PyArray_DATA(py_relative_grid_addresses);
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
  PyArrayObject* py_D_diag;
  PyArrayObject* py_P;
  PyArrayObject* py_Q;
  PyArrayObject* py_A;

  long *D_diag;  /* [3] */
  long (*P)[3];  /* [3][3] */
  long (*Q)[3];  /* [3][3] */
  long (*A)[3];  /* [3][3] */
  long succeeded;

  if (!PyArg_ParseTuple(args, "OOOO",
                        &py_D_diag,
                        &py_P,
                        &py_Q,
                        &py_A)) {
    return NULL;
  }

  D_diag = (long(*))PyArray_DATA(py_D_diag);
  P = (long(*)[3])PyArray_DATA(py_P);
  Q = (long(*)[3])PyArray_DATA(py_Q);
  A = (long(*)[3])PyArray_DATA(py_A);

  succeeded = ksp_get_snf3x3(D_diag, P, Q, A);

  return PyBool_FromLong(succeeded);
}


static PyObject * py_snf_transform_rotations(PyObject *self, PyObject *args)
{
  PyArrayObject* py_transformed_rots;
  PyArrayObject* py_D_diag;
  PyArrayObject* py_Q;
  PyArrayObject* py_rotations;

  long (*transformed_rots)[3][3];
  long *D_diag;  /* [3] */
  long (*Q)[3];  /* [3][3] */
  long (*rotations)[3][3];
  long succeeded, num_rot;

  if (!PyArg_ParseTuple(args, "OOOO",
                        &py_transformed_rots,
                        &py_rotations,
                        &py_D_diag,
                        &py_Q)) {
    return NULL;
  }

  transformed_rots = (long(*)[3][3])PyArray_DATA(py_transformed_rots);
  D_diag = (long(*))PyArray_DATA(py_D_diag);
  Q = (long(*)[3])PyArray_DATA(py_Q);
  rotations = (long(*)[3][3])PyArray_DATA(py_rotations);
  num_rot = PyArray_DIMS(py_rotations)[0];

  succeeded = ksp_snf_transform_rotations(transformed_rots,
                                          rotations, num_rot, D_diag, Q);

  return PyBool_FromLong(succeeded);
}


static PyObject * py_all_grgrid_addresses(PyObject *self, PyObject *args)
{
  PyArrayObject* py_grgrid_address;
  PyArrayObject* py_D_diag;

  long (*grgrid_address)[3];
  long *D_diag;

  if (!PyArg_ParseTuple(args, "OO",
                        &py_grgrid_address,
                        &py_D_diag)) {
    return NULL;
  }

  grgrid_address = (long(*)[3])PyArray_DATA(py_grgrid_address);
  D_diag = (long*)PyArray_DATA(py_D_diag);

  ksp_get_all_grgrid_addresses(grgrid_address, D_diag);

  Py_RETURN_NONE;
}


static PyObject * py_double_grgrid_address(PyObject *self, PyObject *args)
{
  PyArrayObject* py_address_double;
  PyArrayObject* py_address;
  PyArrayObject* py_D_diag;
  PyArrayObject* py_PS;

  long* address_double;
  long* address;
  long* D_diag;
  long* PS;

  if (!PyArg_ParseTuple(args, "OOOO",
                        &py_address_double,
                        &py_address,
                        &py_D_diag,
                        &py_PS)) {
    return NULL;
  }

  address_double = (long*)PyArray_DATA(py_address_double);
  address = (long*)PyArray_DATA(py_address);
  D_diag = (long*)PyArray_DATA(py_D_diag);
  PS = (long*)PyArray_DATA(py_PS);

  ksp_get_double_grgrid_address(address_double,
                                address,
                                D_diag,
                                PS);

  Py_RETURN_NONE;
}


static PyObject * py_grgrid_index(PyObject *self, PyObject *args)
{
  PyArrayObject* py_address;
  PyArrayObject* py_D_diag;

  long* address;
  long* D_diag;
  long grid_index;


  if (!PyArg_ParseTuple(args, "OO",
                        &py_address,
                        &py_D_diag)) {
    return NULL;
  }

  address = (long*)PyArray_DATA(py_address);
  D_diag = (long*)PyArray_DATA(py_D_diag);

  grid_index = ksp_get_grgrid_index(address, D_diag);

  return PyLong_FromLong(grid_index);
}


static PyObject * py_double_grgrid_index(PyObject *self, PyObject *args)
{
  PyArrayObject* py_address_double;
  PyArrayObject* py_D_diag;
  PyArrayObject* py_PS;

  long* address_double;
  long* D_diag;
  long* PS;
  long grid_index;


  if (!PyArg_ParseTuple(args, "OOO",
                        &py_address_double,
                        &py_D_diag,
                        &py_PS)) {
    return NULL;
  }

  address_double = (long*)PyArray_DATA(py_address_double);
  D_diag = (long*)PyArray_DATA(py_D_diag);
  PS = (long*)PyArray_DATA(py_PS);

  grid_index = ksp_get_double_grgrid_index(address_double, D_diag, PS);

  return PyLong_FromLong(grid_index);

}


static PyObject * py_grgrid_address_from_index(PyObject *self, PyObject *args)
{
  PyArrayObject* py_address;
  PyArrayObject* py_D_diag;

  long* address;
  long* D_diag;
  long grid_index;

  if (!PyArg_ParseTuple(args, "OlO",
                        &py_address,
                        &grid_index,
                        &py_D_diag)) {
    return NULL;
  }

  address = (long*)PyArray_DATA(py_address);
  D_diag = (long*)PyArray_DATA(py_D_diag);
  ksp_get_grgrid_address_from_index(address, grid_index, D_diag);

  Py_RETURN_NONE;
}


static PyObject * py_rotate_grgrid_index(PyObject *self, PyObject *args)
{
  PyArrayObject* py_rotation;
  PyArrayObject* py_D_diag;
  PyArrayObject* py_PS;

  long *D_diag;
  long (*rotation)[3];
  long *PS;
  long grid_index;
  long rot_grid_index;


  if (!PyArg_ParseTuple(args, "lOOO",
                        &grid_index,
                        &py_rotation,
                        &py_D_diag,
                        &py_PS)) {
    return NULL;
  }

  rotation = (long(*)[3])PyArray_DATA(py_rotation);
  D_diag = (long*)PyArray_DATA(py_D_diag);
  PS = (long*)PyArray_DATA(py_PS);

  rot_grid_index = ksp_rotate_grgrid_index(grid_index, rotation, D_diag, PS);

  return PyLong_FromLong(rot_grid_index);

}


static PyObject * py_ir_grgrid_map(PyObject *self, PyObject *args)
{
  PyArrayObject* py_ir_grid_indices;
  PyArrayObject* py_rotations;
  PyArrayObject* py_D_diag;
  PyArrayObject* py_PS;

  long *ir_grid_indices;
  long *D_diag;
  long (*rotations)[3][3];
  long *PS;
  long num_rot;

  if (!PyArg_ParseTuple(args, "OOOO",
                        &py_ir_grid_indices,
                        &py_rotations,
                        &py_D_diag,
                        &py_PS)) {
    return NULL;
  }

  ir_grid_indices = (long*)PyArray_DATA(py_ir_grid_indices);
  rotations = (long(*)[3][3])PyArray_DATA(py_rotations);
  num_rot = PyArray_DIMS(py_rotations)[0];
  D_diag = (long*)PyArray_DATA(py_D_diag);
  PS = (long*)PyArray_DATA(py_PS);

  ksp_get_ir_grgrid_map(ir_grid_indices, rotations, num_rot, D_diag, PS);

  Py_RETURN_NONE;
}


static PyObject * py_reciprocal_point_group(PyObject *self, PyObject *args)
{
  PyArrayObject* py_rec_rotations;
  PyArrayObject* py_rotations;
  long is_time_reversal;

  long (*rec_rotations)[3][3];
  long (*rotations)[3][3];
  long num_rot, num_rot_ret;

  if (!PyArg_ParseTuple(args, "OOl",
                        &py_rec_rotations,
                        &py_rotations,
                        &is_time_reversal)) {
    return NULL;
  }

  rec_rotations = (long(*)[3][3])PyArray_DATA(py_rec_rotations);
  rotations = (long(*)[3][3])PyArray_DATA(py_rotations);
  num_rot = PyArray_DIMS(py_rotations)[0];
  num_rot_ret = ksp_get_reciprocal_point_group(rec_rotations,
                                               rotations,
                                               num_rot,
                                               is_time_reversal);
  return PyLong_FromLong(num_rot_ret);
}
