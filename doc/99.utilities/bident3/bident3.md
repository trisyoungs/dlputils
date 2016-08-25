---
title: bident3
brief: Calculate contact types between sites on different molecules, with distance and angle maps
visible: true
taxonomy:
  category: docs
  classification: "Coordination / Contact Numbers"
docroot: /dlputils/docs
template: manpage
---

**bident3** calculates contact types (e.g. monodentate, bidentate etc.) between defined sites on two species, as described for [**bident2**](/dlputils/docs/utilities/bident2), but in addition outputs distance histograms and distance-angle maps for each site and contact type.

## Related Programs
+ To calculate contact distributions, see [**bident2**](/dlputils/docs/utilities/bident2).

## Basic Usage

```
bident3 <HISTORYfile> <OUTPUTfile> <sp1> <sp2> <i> <j> <maxdist> [-frames n] [-discard n] [-site a b] [-dump] [-bin width]
```

Where:
+ `<HISTORYfile>` is the name of the DL_POLY HISTORY file
+ `<OUTPUTfile>` is the name of the DL_POLY OUTPUT file
+ `<sp1>` is the integer index of the species possessing the target 'H' atom sites
+ `<sp2>` is the integer index of the species possessing the 'bidentate' interaction sites
+ `<i>` is the integer index of the first atom of the bidentate binding site on species 2
+ `<j>` is the integer index of the second atom of the bidentate binding site on species 2
+ `<maxdist>` is the maximum allowable distance between sites on the two species

Although this is enough to run **bident3**, at least one 'H' site on species 1 must be defined with the `-site` option.


## Switches

`-bin` _width_

Set the bin width used in distance histograms and the distance component of angle maps to _width_ Angstroms (default = 0.05).

`-discard` _n_

Skip _n_ frames at the beginning of the trajectory before starting the calculation (default = 0)

`-dump`

Print lots of individual contact information to a dump file.

`-frames` _n_

Stop calculating after processing _n_ frames (default = all frames)

`-site` _a_ _b_

Add a site a-b on species 1, where _a_ and _b_ are the local integer indices of the atoms. _a_ represents the 'H'-type atom, and _b_ its bound neighbour. Multiple sites may be defined at once, and will be considered separately in the output.

## Output <a id="output"></a>

Given an input HISTORY file `results.HISu`, then output files are as follows:
`results.bidentNNMM` : all contact information for the defined sites between species NN and MM
`results.single_NNAABB_MMIIJJ.dist` : distance histogram for single contacts between species NN atom AA and species MM atoms II and JJ
`results.bident_NNAABB_MMIIJJ.dist` : distance histogram for bidentate contacts between species NN atom AA and species MM atoms II and JJ
`results.bridge_NNAABB_MMIIJJ.dist` : distance histogram for bridging contacts between species NN atom AA and species MM atoms II and JJ
`results.bifur_NNAABB_MMIIJJ.dist` : distance histogram for bifurcated contacts between species NN atom AA and species MM atoms II and JJ
`results.single_NNAABB_MMIIJJ.surf` : distance-angle map for single contacts between species NN atom AA and species MM atoms II and JJ
`results.bident_NNAABB_MMIIJJ.surf` : distance-angle map for bidentate contacts between species NN atom AA and species MM atoms II and JJ
`results.bridge_NNAABB_MMIIJJ.surf` : distance-angle map for bridging contacts between species NN atom AA and species MM atoms II and JJ
`results.bifur_NNAABB_MMIIJJ.surf` : distance-angle map for bifurcated contacts between species NN atom AA and species MM atoms II and JJ

All `surf` files are three-column format - distance, angle, histogram value - with distinct angles separated by blank lines.
