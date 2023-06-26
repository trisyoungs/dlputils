---
title: cdf
brief: Calculate cylindrical distribution functions
taxonomy:
  category: docs
  classification: "Other Distribution Functions"
docroot: /dlputils/docs
---

**cdf** calculates the distribution function of molecules from a reference vector - as such, it allows one to calculate the distribution function of species confined within cylindrical geometries, for instance, or at the surface of cylindrical objects. Typically, one may use this to calculate distributions of molecules in cylindrical pores, nanotubes etc.

A vector describing the direction along which the cylinder lays is input, along with a point on which this vector lays. For each molecule the minimum distance (i.e. the perpendicular distance) between its centre-of-mass and the defining vector is calculated and binned in histogram form.

Related programs:
+ To calculate cylindrical orientational distribution functions, see [**codf**](/dlputils/docs/utilities/codf).

## Basic Usage

```
cdf <HISTORYfile> <OUTPUTfile> <ox> <oy> <oz> <vx> <vy> <vz>
```

Where:
+ `<HISTORYfile>` is the name of the DL_POLY HISTORY file
+ `<OUTPUTfile>` is the name of the DL_POLY OUTPUT file
+ `<ox> <oy> <oz>` is a point on which the defining vector lays
+ `<vx> <vy> <vz>` is the defining vector of the cylindrical object / pore

## Switches

`-bin` _width_

Specify the histogram bin _width_ to use in the calculation (default = 0.1).

`-compair` _sp_ _i_ _j_

Instead of the centre-of-mass, use the average coordinates of atom indices _i_ and _j_ for species index _sp_. The atom indices are local atom indices for the species in question (i.e. range from 1 to the number of atoms in the molecule).

`-discard` _n_

Discard _n_ frames at the start of the trajectory.

`-frames` _n_

Stop calculating after processing _n_ frames.

`-header` _file_

Read in trajectory header from the specified DL_POLY HIStory _file_. Useful when the desired target frames are present in a restarted (and not appended) HIStory file.

## Normalisation

All CDFs are normalised to reflect the number density of the species, as is standard practice. For systems where the distribution of the species is not homogeneous across the simulation box, the normalisation will be incorrect and must be adjusted by hand.

## Output

Given an input HISTORY file `results.HISu`, then output files are as follows:
`results.cdfNN` : Cylindrical distribution function for species NN about the defining vector


