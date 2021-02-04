# How to use kspclib fortran wrapper

## Fortran wrapper

This fortran wrapper requires c-kspclib library.
https://github.com/atztogo/kspclib/blob/develop/README.md gives
the simple instruction to compile with cmake. Here it is assumed that
the `PREFIX` of the library is located at `../build`.

```
% gfortran -c kspclib_f08.f90
% gfortran example_f08.f90 kspclib_f08.o ../build/lib/libkspc.a -o example
```
