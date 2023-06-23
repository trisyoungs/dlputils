---
title: ardf
brief: Calculate angle-dependent radial distribution functions
visible: true
taxonomy:
  category: docs
  classification: "Radial Distribution Functions"
docroot: /dlputils/docs
template: manpage
---

**ardf** calculates the radial distribution function of sites, binning the results by the angle defined by vectors on each target molecule. The result, then, is a standard [**rdf**](/dlputils/docs/utilities/rdf) split up into functions of individual vector angles. Such functions can be useful when studying, for instance, the variation in approach distance of planar molecules.

Related programs:
+ To calculate centre-of-mass radial distribution functions between species, see [**rdf**](/dlputils/docs/utilities/rdf).
+ To calculate individual atomic partials for one species, see [**rdfaa**](/dlputils/docs/utilities/rdfaa).
+ To calculate site-COM partials, see [**rdfss**](/dlputils/docs/utilities/rdfss).
+ To calculate site-site partials dependent on other sites, see [**rdfdep**](/dlputils/docs/utilities/rdfdep).

## Basic Usage

```
ardf <HISTORYfile> <OUTPUTfile> <sp1> <sp2>
```

Where:
+ `<HISTORYfile>` is the name of the DL_POLY HISTORY file
+ `<OUTPUTfile>` is the name of the DL_POLY OUTPUT file
+ `<sp1>` is the integer index of the species to consider as the central one
+ `<sp2>` is the integer index of the species to compare axes with and form the radial distribution function

Although this is enough to run **ardf**, axis definitions for each of sp1 and sp2 must be provided (see the section on [**Defining Axes**](/dlputils/docs/utilities/pdens#axes) for the [**pdens**](/dlputils/docs/utilities/pdens) code).

## Switches

`-angbin` _width_

Specify the angle histogram bin _width_ to use in the calculation (default = 5.0).

`-axis` _sp_ _x1_ _x2_ _y1_ _y2_

Define the axis system for species index _sp_ from the four supplied atom indices. The atom indices are local atom indices for the species in question (i.e. range from 1 to the number of atoms in the molecule). See the section on [**Defining Axes**](/dlputils/docs/utilities/pdens#axes) for the [**pdens**](/dlputils/docs/utilities/pdens) code for more information.

`-bin` _width_

Specify the histogram bin _width_ to use in the calculation (default = 0.1).

`-cross`

Request calculation and output of cross-axis correlations (i.e. xy, xz, yz etc.).

`-discard` _n_

Discard _n_ frames at the start of the trajectory.

`-frames` _n_

Stop calculating after processing _n_ frames.

`-header` _file_

Read in trajectory header from the specified DL_POLY HIStory _file_. Useful when the desired target frames are present in a restarted (and not appended) HIStory file.

`-restrict` _axis_ _angMin_ _angMax_ _distMin_ _distMax_

By default all surrounding molecules are used to form the ardfs, but the `-restrict` allows a dependency on one axis correlation to be introduced. For instance, one may wish to calculate the distribution of the y and x axis angles if the z axes are aligned more or less parallel. The _axis_ should be the lowercase letter indicating the axis by which to restrict the calculation, which is then followed by the minimum and maximum allowable angles for this like-axis correlation, and those for the distance.


`-sp1site` _atoms_

Defines a site on sp1 that should be used in calculation of distances instead of the origin of the defined axis system. The _atoms_ should be a comma-separated list of local atom indices in the molecule, the avarage position of which will be the new centre.

`-sp2site` _atoms_

Defines a site on sp2 that should be used in calculation of distances instead of the origin of the defined axis system. The _atoms_ should be a comma-separated list of local atom indices in the molecule, the avarage position of which will be the new centre.

## Output <a id="output"></a>

Given an input HISTORY file `results.HISu`, then output files are as follows:
`results.ardfNNMMAABB` : angle RDF between the AA and BB axes of species NN and MM (column 2 = normalised data, column 3 = unnormalised data)
`results.ardfNNMMAABB.sanity` : summed angle RDF between the AA and BB axes of species NN and MM, which should be equivalent to the standard RDF between (written on for AA == BB)

## Normalisation <a id="normalisation"></a>

Since we are splitting up an average over all molecular orientations into ones that depends on the orientation of two vectors, the resulting histogram for a given angle bin must be normalised to account for solid angle effects. Consider two molecules that have their Z axes at 90 degrees to one another. Rotating one vector about the other traces a circle with radius equal to the vector's magnitude. If the angle is now made successively smaller, the circle traced out by rotating the second vector becomes smaller, until at zero degrees (when the vectors overlap) the circle has become an infinitesimally small point in space. This illustrates the 'degeneracy' of the possible ways to orient two free vectors in space - there are more 'ways' of orienting the vectors at higher angles (approaching 90 degrees) than there are at lower angles. For this reason, histograms related to vector angles nearer 90 degrees will have more points contributing to them than will those towards the extremes of zero or 180. As such, a sine function normalisation is applied to the histograms to re-weight them in the correct manner, so that peak heights may be properly compared.

NOTE: For restricted calculations the normalisation is not currently correct!


