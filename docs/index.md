---
title: dlputils User Guide
---

## Overview

**dlputils** is a collection of programs, written in Fortran, which can calculate useful quantitites from a sequence of atomistic configurations. Typically these atomistic configurations represent a molecular dynamics trajectory. These utilities evolved from using the CCP5 [DL_POLY](http://www.scd.stfc.ac.uk/SCD/44516.aspx) to investigate (predominantly) condensed-phase liquid systems.

For visualisation of probability densities calculated with the [`pdens`](pdens) code, consider using [**Aten**](/aten).

## Usage

Most of the codes in **dlputils** package have a similar usage syntax on the command line - at the simplest level, the names of DL_POLY HIStory and OUTput files are the only required arguments. For example:

```
user@pc:~> rdf mysystem.HISu mysystem.OUT
```

Note that the extension of the HIStory file *is* important - passing a default-named HISTORY file makes it difficult to detect whether the file is formatted or unformatted. For this reason, all the **dlputils** assume that a trajectory with the extension _.HISu_ is unformatted, while one with _.HISf_ is formatted. The extension of the OUTput file is not important. It is not necessary to specify the number of molecule types (species) defined in the system, number of atoms in molecules, etc., since this is read from the supplied OUTput file. You can also supply a _0_ (zero) in place of the OUTput file name, in which case you will be asked for the number of molecular species, and the correct number of molecules / atoms for each. The codes are designed to read DL_POLY v2 and DL_CLASSIC HIStory and OUTput files, but making them compatible for DL_POLY v3 (and beyond) should not be a difficult task (modify dlprw as necessary).

Since the first two arguments are mandatory, running any of the codes with no arguments will print out usage information, including a list of available command-line switches which control various relevant parameters in the codes.

In most cases one of more output files are generated, rather than results being written to screen.  The names of these files is always based on the name of the HIStory file, but only if the filename has an extension (it doesn't matter what the extension is). This is removed and replaced by a new, descriptive extension. For example, the **rdf** program produces files called _.rdfNM_, where _N_ and _M_ are the identifying species numbers. So, the above example will produce output files called _mysystem.rdf11_, _mysystem.rdf12_, etc., depending on the number of species in the system.

### Not Using DL_POLY?

Although specific read-write routines for DL_POLY trjectories are the only ones currently implemented, that does not necessarily mean you can't use **dlputils** to analyse your data. As long as you know the unit cell parameters and a standard xyz trajectory, this can be converted to a DL_POLY history file. A suitable OUT file (which is used to describe the frame contents in terms of species, molecules, and atoms) can also be constructed quite easily. See [`xyz2his`](xyz2his) for specific instructions.

## List of Utilities

[`addheader`](addheader) - Add a header to a history file

[`ardf`](ardf) - Calculate angle-dependent radial distribution functions

[`bident2`](bident2) - Calculate contact types between sites on different molecules

[`bident3`](bident3) - Calculate contact types between sites on different molecules, with distance and angle maps

[`cdf`](cdf) - Calculate cylindrical distribution functions

[`cluster`](cluster) - Calculate cluster sizes of molecules

[`clusterab`](clusterab) - Calculate continuous cluster sizes or paths of molecules

[`codf`](codf) - Calculate orientation-dependent cylindrical distribution functions

[`dahist`](dahist) - Calculate distance-angle map between atoms on two species

[`dlp2config`](dlp2config) - Convert a history file into individual CONFIG files

[`dlp2dlp`](dlp2dlp) - Convert between formatted and unformatted HISTORY files

[`dlp2xyzf`](dlp2xyzf) - Convert a history file into a sequential XYZ file (including forces)

[`dlp2xyzs`](dlp2xyzs) - Convert a history file into individual XYZ files

[`dlpfilter`](dlpfilter) - Prune frames from a HIStory file

[`dlpreorder`](dlpreorder) - Reorder / remove species in / from a history file

[`dlpsize`](dlpsize) - Count the frames in a DL_POLY HISTORY file

[`gengg`](gengg) - Generate van der Waals cross terms using geometric mixing rules

[`genlb`](genlb) - Generate van der Waals cross terms using Lorentz-Berthelot mixing rules

[`geom`](geom) - Calculate average geometric parameters

[`getcell`](getcell) - Retrieve and display cell information from history file

[`intertorsion`](intertorsion) - Calculate torsion histogram for a single torsion defined over two species

[`intratorsion`](intratorsion) - Calculate torsion histogram for a single torsion angle defined in one species

[`intratorsion2`](intratorsion2) - Calculate torsion angle map between two torsions defined in one species

[`pdens`](pdens) - Calculate three-dimensional (spatial) probability densities

[`pdensgauss`](pdensgauss) - Apply Gaussian smooth to probability density

[`rdf`](rdf) - Calculate centre-of-mass radial distribution functions

[`rdfaa`](rdfaa) - Calculate atomic partial radial distribution functions involving a single species

[`rdfdep`](rdfdep) - Calculate site-site radial distribution functions dependent on other contacts

[`rdfss`](rdfss) - Calculate partial radial distribution functions from atoms of one species to the centre-of-mass of another

[`sq`](sq) - Calculate neutron-weighted structure factors

[`xyz2his`](xyz2his) - Convert a sequence of xyz trajectory frames to an unformatted DL_POLY HISTORY file

[`zdens`](zdens) - Calculate z-dependent distribution functions

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
