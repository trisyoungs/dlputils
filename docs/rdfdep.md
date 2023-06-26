---
title: rdfdep
brief: Calculate site-site radial distribution functions dependent on other contacts
taxonomy:
  category: docs
  classification: "Radial Distribution Functions"
docroot: /dlputils/docs
---

**rdfdep** allows specific site-site RDFs between two species to be calculated, dependent on the distances of a number of additional specified pairs involving a third species. For instance, in an alcohol/water mixture one may wish to calculate the radial distribution function of those water molecules that are hydrogen bonding directly to the alcohol OH group.

Related programs:
+ To calculate individual atomic partials for one species, see [`rdfaa`](rdfaa).
+ To calculate site-COM partials, see [`rdfss`](rdfss).
+ To calculate centre-of-mass radial distribution functions between species, see [`rdf`](rdf).

## Basic Usage

```
rdfdep <HISTORYfile> <OUTPUTfile> <sp1> <atoms1> <sp2> <atoms2> <depSp>
```

Where:
+ `<HISTORYfile>` is the name of the DL_POLY HISTORY file
+ `<OUTPUTfile>` is the name of the DL_POLY OUTPUT file
+ `<sp1>` is the integer index of the first species
+ `<atoms1>` is a comma-separated list of atoms whose average position form the site of interest on sp1
+ `<sp2>` is the integer index of the second species
+ `<atoms2>` is a comma-separated list of atoms whose average position form the site of interest on sp2
+ `<depSp>` is the integer index of the species for which the 'dependent' contacts with sp1 are tested

## Switches

`-bin` _width_

Specify the histogram bin _width_ to use in the calculation (default = 0.1).

`-dep` _depAtoms_ _sp1atoms_ _minDist_ _maxDist_ _n_

Defines a dependent contact that must be fulfilled in order for the RDF to be calculated. _depAtoms_ and _sp1atoms_ are comma-separated lists of local atom indices on the dependent and first species respectively, defining the sites between which the distance will be tested, with _minDist_ and _maxDist_ defining the acceptable limits. _n_ is the number of times this interaction geometry must be matched for a given sp1 centre to contribute to the final RDF.

`-discard` _n_

Discard _n_ frames at the start of the trajectory.

`-frames` _n_

Stop calculating after processing _n_ frames.

`-header` _file_

Read in trajectory header from the specified DL_POLY HIStory _file_. Useful when the desired target frames are present in a restarted (and not appended) HIStory file.

`-or`

When multiple dependent interactions are specified, require that only one need be fulfilled.

## Output

Given an input HISTORY file `results.HISu`, then output files are as follows:
`results.rdfdepNNMM` : dependent RDF between species NN and MM


