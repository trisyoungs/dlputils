---
title: dlp2config
brief: Convert a history file into individual CONFIG files
taxonomy:
  category: docs
  classification: "History File Conversion / Manipulation"
docroot: /dlputils/docs
---

**dlp2config** converts a DL_POLY history file into one or more CONFIG files

Related programs:
+ To convert a history file into a single sequential XYZ trajectory, see [**dlp2xyzf**](/dlputils/docs/utilities/dlp2xyzf).
+ To extract individual DL_POLY configurations from a history file, see [**dlp2config**](/dlputils/docs/utilities/dlp2config).
+ To convert between formatted and unformatted history files, see [**dlp2dlp**](/dlputils/docs/utilities/dlp2dlp).
+ To extract individual XYZ frames from a history file, see [**dlp2xyzs**](/dlputils/docs/utilities/dlp2xyzs).

## Basic Usage

```
dlp2config <HISTORYfile> <frame>
```

Where:
+ `<HISTORYfile>` is the name of the DL_POLY HISTORY file to be converted
+ `<frame>` is the integer index of the frame to be written, or a negative number specifying the frequency with which to convert multiple frames (e.g. -1 = convert all frames, -10 = convert every 10th frame etc.)

## Switches

None.

