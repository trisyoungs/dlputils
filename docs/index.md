---
title: dlputils User Guide
---

## Overview

**dlputils** is a collection of programs, written in Fortran, which can calculate useful quantitites from a sequence of atomistic configurations. Typically these atomistic configurations represent a molecular dynamics trajectory. These utilities evolved from using the CCP5 [DL_POLY](http://www.scd.stfc.ac.uk/SCD/44516.aspx) to investigate (predominantly) condensed-phase liquid systems.

For visualisation of probability densities calculated with the [**pdens**](/dlputils/docs/utilities/pdens) code, consider using [**Aten**](/aten).

## Compilation

A Fortran90 compiler is required for the bulk of the utilities in **dlputils** - the GNU gfortran compiler works just fine.

**dlputils** uses CMake as its build system - from the top level of the **dlputils** source:

```
bob@pc:~> mkdir build && cd build
```

Then configure:

```
bob@pc:~> cmake ..
```

... and build ...

```
bob@pc:~> cmake --build .
```
