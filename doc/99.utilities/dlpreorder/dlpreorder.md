---
title: dlpreorder
brief: Reorder / remove species in / from a history file
taxonomy:
  category: docs
docroot: /dlputils/docs
template: manpage
---

**dlpreorder** reorders the species in a history file, writing a new history file in the process. If a particular species number is omitted from the list, this species will not appear in the new history file, and so can be used to delete species from a simulation.

## Basic Usage

```
dlpreorder <inputHISTORYfile> <outputHISTORYfile> <species list>
```

Where:
`<inputHISTORYfile>` is the name of the source HISTORY file
`<outputHISTORYfile>` is the target name of the new HIStory file
`<species list>` is a space-separated list of species indices that are to appear in the new file, in the order in which they are given

## Switches

None.

## Output

The new HISTORY file, with the requested name. The new file is always written unformatted.


