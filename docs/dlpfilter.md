---
title: dlpfilter
brief: Prune frames from a HIStory file
taxonomy:
  category: docs
  classification: "History File Conversion / Manipulation"
docroot: /dlputils/docs
---

**dlpfilter** writes a new HIStory file containing a selection of frames from an input HIStory file. The frames to write out can be given either as a text file containing frame indices, one per line, or as a first/skip/last sequence.

## Basic Usage

```
dlpfilter <inputHISTORYfile> <outputHISTORYfile> <framefile | first skip last>
```

Where:
+ `<inputHISTORYfile>` is the name of the source HISTORY file
+ `<outputHISTORYfile>` is the target name of the new HIStory file
`framefile` is the name of a text file containing the desireed frame indices to write in the new file, one per line
`first`, `skip`, `last` are three integer variables stating the starting frame to write, the number of frames to skip inbetween writing frames (use 0 to get all frames), and the final frame to write

## Switches

None.

## Output

The new HISTORY file, with the requested name. The new file is always written unformatted.


