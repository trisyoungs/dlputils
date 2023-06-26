---
title: geom
brief: Calculate average geometric parameters
taxonomy:
  category: docs
  classification: "Intramolecular (Bound) Geometry"
docroot: /dlputils/docs
---

**geom** calculates specific geometric parameters (i.e. bond distances, bond angles, and torsion angles) for a given species, reporting the average, minimum, and maximum values it finds, along with the standard deviation. It can be useful when a quick check on the 'rigidity' of the intramolecular structure is required. Either specific bonds, angles, and torsions may be specified on the command line, or a file containing the required data can be provided.

Related programs:
+ To calculate a single torsion angle histogram, see [**intratorsion**](/dlputils/docs/utilities/intratorsion).
+ To calculate a torsion angle formed between two species, see [**intertorsion**](/dlputils/docs/utilities/intertorsion).
+ To calculate a torsion angle map between two torsions, see [**intratorsion2**](/dlputils/docs/utilities/intratorsion2).
+ To calculate average geometry parameters (not histograms), see [**geom**](/dlputils/docs/utilities/geom).

## Basic Usage

```
geom <HISTORYfile|CONFIGfile> <OUTPUTfile> <sp> <nframes>
```

Where:
+ `<HISTORYfile|CONFIGfile>` is the name of the DL_POLY HISTORY file or, alternatively, a DL_POLY CONFIG file
+ `<OUTPUTfile>` is the name of the DL_POLY OUTPUT file
+ `<sp>` is the integer index of the target species
+ `<nframes>` specifies the number of frames to use in the calculation

## Switches

`-angle` _i_ _j_ _k_

Adds calculation of the geometry of the angle i-j-k.

`-bond` _i_ _j_

Adds calculation of the geometry of the bond i-j.

`-data` _file_

Specifies a _file_ from which bond, angle, and torsion definitions should be read (in addition to any specified on the command line). The format of the file is to contain, one per line, "bond i j", "angle i j k", or "torsion i j k l" where i, j, k, l are atom indices.

`-discard` _n_

Discard _n_ frames at the start of the trajectory.

`-torsion` _i_ _j_ _k_ _l_

Adds calculation of the geometry of the torsion i-j-k-l.

`-warnangle` _maxangle_

Print a warning if any of the target bond angles exceed the supplied _maxangle_.

`-warndist` _maxdist_

Print a warning if any of the target bond distances exceed the supplied _maxdist_.

## Output Files

No output files are produced. Minimum, maximum, average, and standard deviation values are printed to the screen.


