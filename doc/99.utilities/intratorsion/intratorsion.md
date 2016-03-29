---
title: intratorsion
brief: Calculate torsion histogram for a single torsion angle defined in one species
visible: true
taxonomy:
  category: docs
  classification: "Intramolecular (Bound) Geometry"
docroot: /dlputils/docs
template: manpage
---

**intratorsion** calculates a histogram of angles, from -180 to +180, for a single torsion angle in a given species.

Related programs:
+ To calculate a torsion angle formed between two species, see [**intertorsion**](/dlputils/docs/intertorsion).
+ To calculate a torsion angle map between two torsions, see [**intratorsion2**](/dlputils/docs/intratorsion2).
+ To calculate average geometry parameters (not histograms), see [**geom**](/dlputils/docs/geom).

## Basic Usage

```
intratorsion <HISTORYfile> <OUTPUTfile> <sp> <i> <j> <k> <l>
```

Where:
`<HISTORYfile>` is the name of the DL_POLY HISTORY file
`<OUTPUTfile>` is the name of the DL_POLY OUTPUT file
`<sp>` is the integer index of the target species
`<i>` is the integer index of the atom 'i' in the torsion i-j-k-l
`<j>` is the integer index of the atom 'j' in the torsion i-j-k-l
`<k>` is the integer index of the atom 'k' in the torsion i-j-k-l
`<l>` is the integer index of the atom 'l' in the torsion i-j-k-l

## Switches

None.

## Output Files

The torsion angle histogram is written to a file named after the input `<HISTORYfile>`, ending with `<i>-<j>-<k>-<l>.tors`.


