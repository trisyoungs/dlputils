---
title: intertorsion
brief: Calculate torsion histogram for a single torsion defined over two species
visible: true
taxonomy:
  category: docs
  classification: "Intramolecular (Bound) Geometry"
docroot: /dlputils/docs
template: manpage
---

**intertorsion** calculates a histogram of torsion angles, from -180 to +180, for a single torsion angle defined by two atoms on one species and two atoms on another (not necessarily different) species. As such, it provise the possibility to analyse torsion angles formed by, for instance, hydrogen bond donor and acceptor groups. As well as the full torsion histogram, histograms for the two angles and the inter-species distance are also provided, the latter in the form of a 2D map.

Related programs:
+ To calculate a single torsion angle histogram, see [**intratorsion**](/dlputils/docs/utilities/intratorsion).
+ To calculate a torsion angle map between two torsions, see [**intratorsion2**](/dlputils/docs/utilities/intratorsion2).
+ To calculate average geometry parameters (not histograms), see [**geom**](/dlputils/docs/utilities/geom).

## Basic Usage

```
intratorsion <HISTORYfile> <OUTPUTfile> <sp1> <i> <j> <sp2> <k> <l> [maxjk]
```

Where:
+ `<HISTORYfile>` is the name of the DL_POLY HISTORY file
+ `<OUTPUTfile>` is the name of the DL_POLY OUTPUT file
+ `<sp1>` is the integer index of the target species containing atoms i and j
+ `<i>` is the integer index of the atom 'i', in species 1, in the torsion i-j-k-l
+ `<j>` is the integer index of the atom 'j', in species 1, in the torsion i-j-k-l
+ `<sp2>` is the integer index of the target species containing atoms k and l
+ `<k>` is the integer index of the atom 'k', in species 2, in the torsion i-j-k-l
+ `<l>` is the integer index of the atom 'l', in species 2, in the torsion i-j-k-l
...and optionally...
[`maxjk`] is the maximum allowable distance between atoms j and k (if not given, any distance is allowed).

## Switches

None.

## Output Files

The main torsion angle histogram is written to a file named after the input + `<HISTORYfile>`, ending with `<i>-<j>-<k>-<l>.ijkl`. Other output files are:
+ `<HISTORYfile>.<j>-<k>.jk`, containing the jk distance histogram
+ `<HISTORYfile>.<i>-<j>-<k>.ijk`, containing the ijk angle histogram
+ `<HISTORYfile>.<j>-<k>-<l>.jkl`, containing the jkl angle histogram
+ `<HISTORYfile>.<i>-<j>-<k>-<l>.ijkl_jk`, containing the jk / ijkl (distance / torsion) 2D histogram map


