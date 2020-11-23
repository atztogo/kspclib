#include <Python.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include <kspclib.h>

#if (PY_MAJOR_VERSION < 3) && (PY_MINOR_VERSION < 6)
#define PYUNICODE_FROMSTRING PyString_FromString
#else
#define PYUNICODE_FROMSTRING PyUnicode_FromString
#endif

static PyObject * py_kspclib(PyObject *self, PyObject *args);

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
  {"kspclib", py_kspclib, METH_VARARGS, "Kspclib"},
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

static PyObject * py_kspclib(PyObject *self, PyObject *args)
{
  PyArrayObject* py_A;
  PyArrayObject* py_P;
  PyArrayObject* py_Q;

  int (*A)[3];
  int (*P)[3];
  int (*Q)[3];

  A = NULL;
  P = NULL;
  Q = NULL;

  if (!PyArg_ParseTuple(args, "OOO", &py_A, &py_P, &py_Q)) {
    return NULL;
  }

  A = (int(*)[3])PyArray_DATA(py_A);
  P = (int(*)[3])PyArray_DATA(py_P);
  Q = (int(*)[3])PyArray_DATA(py_Q);

  kspclib(A, P, Q);

  Py_RETURN_NONE;
}
