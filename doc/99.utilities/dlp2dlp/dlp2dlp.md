---
title: dlp2dlp
brief: Convert between formatted and unformatted HISTORY files
visible: true
taxonomy:
  category: docs
  classification: "History File Conversion / Manipulation"
docroot: /dlputils/docs
template: manpage
---

**dlp2dlp** allows conversion of a formatted HISTORY file to its unformatted analogue, and vice versa. Note that converting unformatted HISTORY files created on an architecture different from the one you are running **dlp2dlp** will fail.

## Basic Usage

```
dlp2dlp <HISTORYfile> <format> [nframes]
```

Where:
+ `<HISTORYfile>` is the target name of the DL_POLY HISTORY file to be written
+ `<format>` is the desired format of the output file, _0_ for unformatted, or _1_ for formatted
`[nframes]` (optional) is the number of frames to convert (default = all)

## Switches

None.

## Output

The new HISTORY file is called `converted.HISf` or `converted.HISu`, depending on the format requested.


