# kspclib

C library for handling regular grid sampling.

## How to compile with cmake

```
% mkdir _build
% cd _build
% cmake ..
% make
% make install (probably installed under /usr/local)
```

Or to install under /some/where

```
% mkdir _build
% cd _build
% cmake -DCMAKE_INSTALL_PREFIX="" ..
% make
% make DESTDIR=/some/where install
```
