---
title: clusterab
brief: Calculate continuous cluster sizes or paths of molecules
taxonomy:
  category: docs
docroot: /dlputils/docs
header_class: alt
---

**clusterab** calculates the sizes of clusters (or semi-continuous networks) of molecules in a similar manner to [a]**cluster**,dlputils,manual,cluster[/a], but focuses specifically on 'A-B' contacts between the various molecules. Thus, rather than giving the cluster size of a group of molecules based on a single common site between them, **clusterab** allow true networks to be probed. The best example of this is a hydrogen bonding network, where the hydrogen (A) atoms of one molecule donate a hydrogen bond to an oxygen or nitrogen (B) atom on another (possibly different) molecule. Also, multiple A and B sites may be specified for each molecule type involved, allowing diols, for instance, to be treated in the same manner.

Since the present calculated quantity may involve more than one species, the resulting cluster sizes the sum of the involved molecules of both types. One exception to this is the calculation of a cluster size when then surrounding molecules of interest do not include the central species type - in this case, cluster sizes of zero are possible, since the central molecule may not be interacting with any of the other molecule types.

The normal action of **clusterab** is to allow paths (i.e. A-B interactions) involving as many molecules as possible, and to allow each molecule in the system to belong to exactly one cluster/path. Sometimes this is not the desired quantity - if, for instance, one was interested in the average number of a molecule type that were clustering around an AB site of a (different) central species molecule, then each molecule may be associated (in)directly to more than one species 1 AB site. In these situations, the `-individual` option may be provided, which has the effect of 'forgetting' the previous path before the next molecule of the central species is considered in the calculation. The allowable path length in either case can be restricted with the `-maxpath` switch, allowing 'extended' coordination shells to be investigated.

## Basic Usage

```
clusterab <HISTORYfile> <OUTPUTfile> <sp> <othersp> <maxdist>
```

Where:
`<HISTORYfile>` is the name of the DL_POLY HISTORY file
`<OUTPUTfile>` is the name of the DL_POLY OUTPUT file
`<sp>` is the index of the central species, from which the cluster calculation begins
`<othersp>` is a comma-separated list of indices to test use in the calculation of cluster sizes
`<maxdist>` is the maximum allowable distance between any A-B site pair for them to be considered 'interacting'.

## Switches

`-ab` _sp_ _aAtoms_ _bAtoms_
Specify A and B sites for the species index _sp_. _aAtoms_ and _bAtoms_ are comma-separated lists of local atom indices on the species. The donor or acceptor sites can be assigned to either A or B sites, but you must be consistent across all species.

`-discard` _n_
Discard _n_ frames at the start of the trajectory.

`-frames` _n_
Stop calculating after processing _n_ frames.

`-header` _file_
Read in trajectory header from the specified DL_POLY HIStory _file_. Useful when the desired target frames are present in a restarted (and not appended) HIStory file.

`-individual`
Treat each central species molecule as a separate calculation, allowing any of the other species to be associated to its cluster.

`-maxpath` _n_
Sets the maximum number of sequential A-B interactions to allow before stopping the recursive search. A value of zero will allow no clusters to be formed, a value of 1 will restrict to 'immediate' neighbours, etc. The default is to allow paths of any length.

## Output

Given an input HISTORY file **results.HISu**, then output files are as follows:
**results.clusterabSP** : cluster information for species SP

</page>
