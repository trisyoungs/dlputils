---
title: codf
brief: Calculate orientation-dependent cylindrical distribution functions
taxonomy:
  category: docs
  classification: "Other Distribution Functions"
docroot: /dlputils/docs
---

**codf** acts in a similar manner to [**cdf**](/dlputils/docs/utilities/cdf), calculating cylindrical distribution functions of molecules. However, instead of the normalised 2D equivalent to the [**rdf**](/dlputils/docs/utilities/rdf), **codf** gives information on the orientation of species as a function of distance from the defining vector, given defining sets of axes on the individual species.

A vector describing the direction along which the cylinder lays is input, along with a point on which this vector lays. For each molecule the minimum distance (i.e. the perpendicular distance) between its centre-of-mass and the defining vector is calculated and the dot product between the perpendicular vector and each axis of the molecule is summed into distance bins. The resulting output is the average of the dot product of each molecular axis as a function of distance from the defining vector, along with the 2D maps of angle vs distance for each of the three angles formed between the perpendicular vector and the molecule axes.

Related programs:
+ To calculate standard cylindrical distribution functions, see [**cdf**](/dlputils/docs/utilities/cdf).

## Basic Usage

```
codf <HISTORYfile> <OUTPUTfile> <ox> <oy> <oz> <vx> <vy> <vz> <targetsp>
```

Where:
+ `<HISTORYfile>` is the name of the DL_POLY HISTORY file
+ `<OUTPUTfile>` is the name of the DL_POLY OUTPUT file
+ `<ox> <oy> <oz>` is a point on which the defining vector lays
+ `<vx> <vy> <vz>` is the defining vector of the cylindrical object / pore
+ `<targetsp>` is a comma-separated list of indices specifying which species to consider

## Switches

`-anglebin` _width_

Specify the angle histogram bin _width_ to use in the calculation (default = 1.0).

`-axis` _sp_ _x1_ _x2_ _y1_ _y2_

Define the axis system for species index _sp_ from the four supplied atom indices. Axes must be defined for each + `<targetsp>` listed. See the [**pdens** manual](/dlputils/docs/utilities/pdens#axes) for more information on defining axes.

`-bin` _width_

Specify the distance histogram bin _width_ to use in the calculation (default = 0.1).

`-compair` _sp_ _i_ _j_

Instead of the centre-of-mass, use the average coordinates of atom indices _i_ and _j_ for species index _sp_. The atom indices are local atom indices for the species in question (i.e. range from 1 to the number of atoms in the molecule).

`-discard` _n_

Discard _n_ frames at the start of the trajectory.

`-frames` _n_

Stop calculating after processing _n_ frames.

`-header` _file_

Read in trajectory header from the specified DL_POLY HIStory _file_. Useful when the desired target frames are present in a restarted (and not appended) HIStory file.

## Output

Given an input HISTORY file `results.HISu`, then output files are as follows:
`results.codfNN` : Cylindrical orientation distribution function for species NN about the defining vector (columns 2-4 are standard averages, columns 5-7 are averages of the absolute dot products, column 8 is the number of molecules used in the average)
`results.codfNNAA` : Distance / angle map for axis AA of species NN


