---
title: bident2
brief: Calculate contact types between sites on different molecules
taxonomy:
  category: docs
  classification: "Coordination / Contact Numbers"
docroot: /dlputils/docs
---

**bident2** calculates contact types (e.g. monodentate, bidentate etc.) between defined sites on two species. It came about initially because of the need to analyse simulations comprising acetate anions (hence the 'bident'-ate name of the program) where the type of binding between the acetate and acidis hydrogen sites on a second species was unknown. Because of this the code is quite specific, but can be 'abused' to analyse contact numbers for systems where the primary interaction site is not potentially bidentate.

The code analyses all contacts less than a specified distance, grouping them into monodentate (where an 'H' site is within range of one of the two sites of the second species), bidentate (where an 'H' site is within range of both sites of the second species), bridging (where two 'H' sites are bridged by the second species, one contact for each of its sites), bifurcated (where one 'H' site has a contact with both sites of the second species simultaneously), and multiple (which encompasses all other cases that are not so easy to interpret or name).

The initial intention of **bident2** was to analyse hydrogen bonding contacts between acidic H sites of a cation and the bidentate oxygen moeity of an anion, but as already stated the code can be lightly 'abused'. Firstly, there is  no need for the two sites on the second species to be situated such that, together, they might be considered a potentially 'bidentate' binding site. The program will run just fine, but of course no bidentate (or likely any other multi-contact interactions) will be found, but the single contacts will be caught. In this way **bident2** can be used as a probe for contact numbers between individual binding sites.

## Related Programs
+ To calculate contact distributions and output distance / angle maps, see [**bident3**](/dlputils/docs/utilities/bident3).

## Basic Usage

```
bident2 <HISTORYfile> <OUTPUTfile> <sp1> <sp2> <i> <j> <maxdist> [-frames n] [-discard n] [-site a b] [-dump]
```

Where:
+ `<HISTORYfile>` is the name of the DL_POLY HISTORY file
+ `<OUTPUTfile>` is the name of the DL_POLY OUTPUT file
+ `<sp1>` is the integer index of the species possessing the target 'H' atom sites
+ `<sp2>` is the integer index of the species possessing the 'bidentate' interaction sites
+ `<i>` is the integer index of the first atom of the bidentate binding site on species 2
+ `<j>` is the integer index of the second atom of the bidentate binding site on species 2
+ `<maxdist>` is the maximum allowable distance between sites on the two species

Although this is enough to run **bident2**, at least one 'H' site on species 1 must be defined with the `-site` option.


## Switches

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

This outout file lists each defined `-site` in its own column, detailing the total number of contacts found (regardless of their style) for each, followed by a breakdown of the individual contact styles by number and percentages of group and total contact numbers.

