---
title: cluster
brief: Calculate cluster sizes of molecules
visible: true
taxonomy:
  category: docs
  classification: "Coordination / Contact Numbers"
docroot: /dlputils/docs
template: manpage
---

**cluster** calculates the sizes of clusters (or semi-continuous networks) of molecules, given some distance criteria and a specific site on the molecule. By default, the centre-of-mass is used as the point of interest, but this can be defined to be any point on the molecule. Typically one may choose, for instance, the oxygen atom of an OH group as the site of interest, and set the cutoff to be sensible with respect to O-O distances in typical hydrogen bonds (e.g. 3.1 Angstroms or slightly less).

## Basic Usage

```
cluster <HISTORYfile> <OUTPUTfile> <sp> <maxdist>
```

Where:
+ `<HISTORYfile>` is the name of the DL_POLY HISTORY file
+ `<OUTPUTfile>` is the name of the DL_POLY OUTPUT file
+ `<sp>` is the index of the target species
+ `<maxdist>` is the maximum allowable distance between the COM / point on the molecule

## Switches

`-com` _indices_

Instead of the centre-of-mass, use the average coordinates of the (comma-separated) list of atom _indices_ provided as the site of interest on the species. The atom indices are local atom indices for the species in question (i.e. range from 1 to the number of atoms in the molecule).

`-discard` _n_

Discard _n_ frames at the start of the trajectory.

`-frames` _n_

Stop calculating after processing _n_ frames.

`-header` _file_

Read in trajectory header from the specified DL_POLY HIStory _file_. Useful when the desired target frames are present in a restarted (and not appended) HIStory file.

## Output

Given an input HISTORY file `results.HISu`, then output files are as follows:
`results.clusterSP` : cluster information for species SP

