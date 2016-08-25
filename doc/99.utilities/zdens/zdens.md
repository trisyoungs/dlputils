---
title: zdens
brief: Calculate z-dependent distribution functions
visible: true
taxonomy:
  category: docs
  classification: "Other Distribution Functions"
docroot: /dlputils/docs
template: manpage
---

**zdens** calculates the distribution of the centres-of-mass (averaged over species) of all molecules in the system as a function of z-coordinate.

**NOTE**: zdens will not work correctly if the axes are not orthogonal.

Related programs:
+ To calculate radial distribution functions, see [**rdf**](/dlputils/docs/utilities/rdf).

## Basic Usage

```
zdens <HISTORYfile> <OUTPUTfile> <base Z>
```

Where:
+ `<HISTORYfile>` is the name of the DL_POLY HISTORY file
+ `<OUTPUTfile>` is the name of the DL_POLY OUTPUT file
+ `<base Z>` is the z-coordinate to use as the 'zero' point for the distribution (e.g. it could be the location of a surface, or just zero to get the distribution relative to the box itself)

## Switches

`-bin` _width_

Specify the histogram bin _width_ to use in the calculation (default = 0.1).

`-comatom` _sp_ _i_

Instead of the centre-of-mass, use the specifed atom index _i_ as the reference point for species index _sp_. The atom indices are local atom indices for the species in question (i.e. range from 1 to the number of atoms in the molecule).

`-discard` _n_

Discard _n_ frames at the start of the trajectory.

`-frames` _n_

Stop calculating after processing _n_ frames.

`-header` _file_

Read in trajectory header from the specified DL_POLY HIStory _file_. Useful when the desired target frames are present in a restarted (and not appended) HIStory file.

`-symm` _centreZ_

Symmetrise the resulting distribution about the supplied _centreZ_ coordinate.

## Normalisation

Distribution functions are normalised to frame count alone.

## Output

Given an input HISTORY file `results.HISu`, then output files are as follows:
`results.zdnNN` : z-density distribution for species NN


