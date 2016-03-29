---
title: intratorsion2
brief: Calculate torsion angle map between two torsions defined in one species
visible: true
taxonomy:
  category: docs
  classification: "Intramolecular (Bound) Geometry"
docroot: /dlputils/docs
template: manpage
---

Calculates a 2D map (histogram) between two torsion angles on the same species.

Related programs:
+ To calculate a torsion angle formed between two species, see [**intertorsion**](/dlputils/docs/intertorsion).
+ To calculate a single torsion angle histogram, see [**intratorsion**](/dlputils/docs/intratorsion).
+ To calculate average geometry parameters (not histograms), see [**geom**](/dlputils/docs/geom).

## Basic Usage

```
intratorsion2 <HISTORYfile> <OUTPUTfile> <sp> <i1> <j1> <k1> <l1> <i2> <j2> <k2> <l2>
```

Where:
`<HISTORYfile>` is the name of the DL_POLY HISTORY file
`<OUTPUTfile>` is the name of the DL_POLY OUTPUT file
`<sp>` is the integer index of the target species
`<i1>` is the integer index of the atom 'i' in the first torsion i-j-k-l
`<j1>` is the integer index of the atom 'j' in the first torsion i-j-k-l
`<k1>` is the integer index of the atom 'k' in the first torsion i-j-k-l
`<l1>` is the integer index of the atom 'l' in the first torsion i-j-k-l
`<i2>` is the integer index of the atom 'i' in the second torsion i-j-k-l
`<j2>` is the integer index of the atom 'j' in the second torsion i-j-k-l
`<k2>` is the integer index of the atom 'k' in the second torsion i-j-k-l
`<l2>` is the integer index of the atom 'l' in the second torsion i-j-k-l

## Switches

None.

## Output Files

The torsion angle histogram is written to a file named after the input `<HISTORYfile>`, ending with `<i1>-<j1>-<k1>-<l1>_<i2>-<j2>-<k2>-<l2>.tors2`.

