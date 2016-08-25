---
title: dahist
brief: Calculate distance-angle map between atoms on two species
visible: true
taxonomy:
  category: docs
docroot: /dlputils/docs
template: manpage
---

**dahist** calculates a distance-angle map for an A...B-C interaction, giving both the A-B-C angle map and the A-B distance histrogram as output.

Related programs:
+ To calculate individual atomic partials for one species, see [**rdfaa**](/dlputils/docs/utilities/rdfaa).

## Basic Usage

```
dahist <HISTORYfile> <OUTPUTfile> <sp1> <sp2> <maxdist> [-triplet a b c] [-frames n] [-mindist d]
```

Where:
+ `<HISTORYfile>` is the name of the DL_POLY HISTORY file
+ `<OUTPUTfile>` is the name of the DL_POLY OUTPUT file
+ `<sp1>` is the integer index of the species possessing atom 'a' in the interaction
+ `<sp2>` is the integer index of the species possessing atoms 'b' and 'c' in the interaction
+ `<maxdist>` is the maximum allowable distance for a-b for it to be considered

Although this is enough to run **dahist**, at least one triplet of atoms muse be specified with the `-triplet` switch.


## Switches

`-frames` _n_

Stop calculating after processing _n_ frames.

`-mindist` _d_

Specifies a minimum allowable distance for the atoms a-b, below which they will not be considered (default = 0.0).

`-triplet` _a_ _b_ _c_

Add a triplet a-b-c to the calculation. Note that if multiple triplets are defined these contribute to the same distance-angle map - if separate maps are required for each triplet the code must be called individually for each.

## Output <a id="output"></a>

Given an input HISTORY file `results.HISu`, then output files are as follows:
`results.dahistNNAAMMBBCC.hist` : histogram of the a-b distances found in the calculation
`results.dahistNNAAMMBBCC.surf` : distance angle map of the a-b-c interaction in three column format (distance, angle, histogram value) with a blank line separating each angle

