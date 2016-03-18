---
title: genlj
brief: Generate van der Waals cross terms using Lorentz-Berthelot mixing rules
taxonomy:
  category: docs
docroot: /dlputils/docs
template: manpage
---

Given a file containing three columns (atom type names, epsilon and sigma values, separated by spaces - essentially the format used by the DL_POLY FIELD file in the 'vdw' section) **genlj** will output the full list of self and cross-terms, using geometric mixing for epsilon and arithmetic mixing for sigma.

## Basic Usage

```
genlj <PARAMfile>
```

Where:
`<PARAMfile>` is the name of the file containing atomtype name, epsilon, and sigma values, separated by spaces

## Switches

None.

## Output

The full list of cross terms is written to stdout.

