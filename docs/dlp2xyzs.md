---
title: dlp2xyzs
brief: Convert a history file into individual XYZ files
taxonomy:
  category: docs
  classification: "History File Conversion / Manipulation"
docroot: /dlputils/docs
---

**dlp2xyzs** converts a DL_POLY history file into a separate XYZ files.

Related programs:
+ To convert a history file into a single sequential XYZ trajectory, see [`dlp2xyzf`](dlp2xyzf).
+ To extract individual DL_POLY configurations from a history file, see [`dlp2config`](dlp2config).
+ To convert between formatted and unformatted history files, see [`dlp2dlp`](dlp2dlp).

## Basic Usage

```
dlp2xyzf <HISTORYfile> <interval>
```

Where:
+ `<HISTORYfile>` is the name of the DL_POLY HISTORY file to be converted
+ `<interval>` is the frequency with which to output XYZ frames of the trajectory (1 = write every frame, 10 = write every 10th frame etc.)

## Switches

None.

