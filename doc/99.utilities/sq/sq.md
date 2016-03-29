---
title: sq
brief: Calculate neutron-weighted structure factors
visible: true
taxonomy:
  category: docs
docroot: /dlputils/docs
template: manpage
---

**sq** calculates the neutron-weighted F(Q) (coherent bound cross-section I(Q), differential cross section, static structure factor etc.) of the system, using partials in line with the Faber-Ziman formulism. 

## Arguments

TODO

Basic usage:

```
sq <HISTORYfile> <OUTPUTfile>
```

Where:
`<HISTORYfile>` is the name of the DL_POLY HISTORY file
`<OUTPUTfile>` is the name of the DL_POLY OUTPUT file

Switches

TODO

## Atomtype / Isotope Definitions

**sq** requires that a file is provided which contains atom type name to element/isotope mappings. By default, this is expected to be called _lengths.dat_, but an alternative name can be specified with the `-lengths` switch. This file is needed to enable the code to unambiguously distinguish which atoms correspond to which elements (something which is not always obvious from the atom type names used in forcefields). The file also allows atom type names to be mapped to a new, alternative name, or to combine several atom type names into one single type name in the program's output. Currently, only a basic set of elements and their neutron scattering lengths are defined in the actual code, meaning that exotic (read: 'interesting') systems may fail. However, adding the relevant data to the code should be a relatively simple task. Until the entire Sears table is added, the onus is on the end user to perform this task.

The first line of the file should consist of a single integer specifying the number of lines of data to expect in the remainder of the file. After this, each line should contain three items, separated by tabs or spaces; the original atom type name (as used in the HISTORY file), the new name to map to, and finally the element to be used for this atom type (or D for deuterium). For example:

```
5
CA	C	C
CB	C	C
HA	HA	H
HB	HB	D
N	N	N
```

This example states that both carbon types should be considered as equivalent (and just be named 'C') while all other aom types retain their original names. Note, however, that _HA_ is treated as hydrogen and _HB_ is treated as deuterium. The element symbol may also be _X_ which indicates a scattering length of zero should be used for this atom type.  The mappings are important since the number and names of the (optional) partial S(Q) output files (requested with the `-partials` switch) depends on the number and names of 'remapped' names (second column) in the file. All names are case-sensitive.


