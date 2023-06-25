---
title: gengg
brief: Generate van der Waals cross terms using geometric mixing rules
visible: true
taxonomy:
  category: docs
  classification: "Miscellaneous"
docroot: /dlputils/docs
template: manpage
---

Given a file containing three columns (atom type names, epsilon and sigma values, separated by spaces - essentially the format used by the DL_POLY FIELD file in the 'vdw' section) **gengg** will output the full list of self and cross-terms, using geometric mixing rules for the parameters.

## Basic Usage

```
gengg <PARAMfile>
```

Where:
+ `<PARAMfile>` is the name of the file containing atomtype name, epsilon, and sigma values, separated by spaces

## Switches

None.

## Output

The full list of cross terms is written to stdout.


