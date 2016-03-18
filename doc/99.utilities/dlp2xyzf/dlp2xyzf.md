---
title: dlp2xyzf
brief: Convert a history file into a sequential XYZ file (including forces)
taxonomy:
  category: docs
docroot: /dlputils/docs
template: manpage
---

**dlp2xyzf** converts a DL_POLY history file into a sequential XYZ file, with atomic positions and forces contained in the file.

Related programs:
+ To extract individual DL_POLY configurations from a history file, see [**dlp2config**](/dlputils/docs/dlp2config).
+ To convert between formatted and unformatted history files, see [**dlp2dlp**](/dlputils/docs/dlp2dlp).
+ To extract individual XYZ frames from a history file, see [**dlp2xyzs**](/dlputils/docs/dlp2xyzs).

## Basic Usage

```
dlp2xyzf <HISTORYfile> <nframes | 'all'>
```

Where:
`<HISTORYfile>` is the name of the DL_POLY HISTORY file to be converted
`<nframes | 'all'>` is either the integer number of frames to convert, or the word 'all' to convert all frames in the trajectory.

## Switches

None.

