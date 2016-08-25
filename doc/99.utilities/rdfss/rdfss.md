---
title: rdfss
brief: Calculate partial radial distribution functions from atoms of one species to the centre-of-mass of another
visible: true
taxonomy:
  category: docs
  classification: "Radial Distribution Functions"
docroot: /dlputils/docs
template: manpage
---

**rdfss** calculates calculates atomic, partial radial distribution functions between all individual atoms of one species, to the centre-of-mass (or another specified point) of another).

Related programs:
+ To calculate centre-of-mass radial distribution functions between species, see [**rdf**](/dlputils/docs/utilities/rdf).
+ To calculate individual atomic partials for one species, see [**rdfaa**](/dlputils/docs/utilities/rdfaa).
+ To calculate site-site partials dependent on other sites, see [**rdfdep**](/dlputils/docs/utilities/rdfdep).

## Basic Usage

```
rdfss <HISTORYfile> <OUTPUTfile> <sp1> <sp2>
```

Where:
+ `<HISTORYfile>` is the name of the DL_POLY HISTORY file
+ `<OUTPUTfile>` is the name of the DL_POLY OUTPUT file
+ `<sp1>` is the integer index of the species for which all atoms are considered
+ `<sp2>` is the integer index of the species for which the centre-of-mass is considered

## Switches

`-bin` _width_

Specify the histogram bin _width_ to use in the calculation (default = 0.1).

`-compair` _i_ _j_

Instead of the centre-of-mass, use the average coordinates of atom indices _i_ and _j_ for the second species. The atom indices are local atom indices for the species in question (i.e. range from 1 to the number of atoms in the molecule).

`-discard` _n_

Discard _n_ frames at the start of the trajectory.

## Output

Given an input HISTORY file `results.HISu`, then output files are as follows:
`results.ssrdfCC_SP_IINN` : partial RDF between atom II (name NN) of species 1 (SP) to species 2 (CC)


