---
title: Compilation
brief: Notes and instructions for  compiling dlputils
taxonomy:
  category: docs
docroot: /dlputils/docs
template: manpage
---

A Fortran90 compiler is required for the bulk of the utilities in **dlputils** - the GNU gfortran compiler works just fine. **dlputils** uses the standard configure / make toolchain or CMake, so to initialise the build process simply perform your choice the following steps (if in doubt, go with the configure / make option).

## Autotools (Configure / Make)

If you retrieved the source code via svn or git, you will need to run autogen.sh first (note: this requires the **autotools** and **automake** packages to be installed on your system):

```
bob@pc:~> ./autogen.sh
```

Now, run the main configure script:

```
bob@pc:~> ./configure
```

...followed by...

```
bob@pc:~> make
```

This is enough to compile most of the codes. There are a few utilities (**acf**, **dacf**, **sq**, and **totalgr**) which are computationally quite intensive and require MPI (for instance [OpenMPI](http://www.open-mpi.org)), and which must be compiled separately. Provided a suitable mpif90 executable is somewhere in your path, compiling these individually is as simple as:

```
bob@pc:~> make sq
```

## CMake

An out-of-tree build is recommended. From an empty build directory of your choice (which could simply be a new subdirectory within the dlputils source tree), run:

```
bob@pc:~> cmake /path/to/dlputils
```

For example:

```
bob@pc:~> cmake ../
```

...if you created a new directory in the dlputils source tree (so the Fortran sources are one directory up from your current location).

On Windows you will need to specify the generator to use - for the MinGW toolchain, use:

```
bob@pc:~> cmake -G "MinGW Makefiles" ../
```


