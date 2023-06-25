---
title: rdfaa
brief: Calculate atomic partial radial distribution functions involving a single species
visible: true
taxonomy:
  category: docs
  classification: "Radial Distribution Functions"
docroot: /dlputils/docs
template: manpage
---

**rdfaa** calculates calculates atomic, partial radial distribution functions between specific atoms of a single species. Contributions from atom pairs in the same species are omitted, unless specifically requested.

Related programs:
+ To calculate centre-of-mass radial distribution functions between species, see [**rdf**](/dlputils/docs/utilities/rdf).
+ To calculate site-COM partials, see [**rdfss**](/dlputils/docs/utilities/rdfss).
+ To calculate site-site partials dependent on other sites, see [**rdfdep**](/dlputils/docs/utilities/rdfdep).

## Basic Usage

```
rdfaa <HISTORYfile> <OUTPUTfile> <sp>
```

Where:
+ `<HISTORYfile>` is the name of the DL_POLY HISTORY file
+ `<OUTPUTfile>` is the name of the DL_POLY OUTPUT file
+ `<sp>` is the integer index of the target species

## Switches

`-discard` _n_

Discard _n_ frames at the start of the trajectory.

`-intra`

Intramolecular (same molecule) contributions to the partial RDFs will be included in the final results.

`-frames` _n_

Stop calculating after processing _n_ frames.

`-pair` _i_ _j_

Add the atom indices i-j to the list of partial RDFs to calculate for the specified species. The atom indices are local atom indices for the species in question (i.e. range from 1 to the number of atoms in the molecule).

## Output

Given an input HISTORY file `results.HISu`, then output files are as follows:
`results.rdfaaSP_II_JJ` : partial RDF between atom II and JJ in species SP


