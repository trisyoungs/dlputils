---
title: addheader
brief: Add a header to a history file
taxonomy:
  category: docs
docroot: /dlputils/docs
template: manpage
---

For a couple of reasons, trajectory files generated from DL_POLY occasionally have no header - **addheader** allows one to be copied from an existing history file, and prepended to the target history file.

## Basic Usage

```
addheader <HEADERfile> <input HISTORYfile> <output HISTORYfile>
```

Where:
`<HEADERfile>` is the name of the DL_POLY HISTORY file from which the header is to be copied
`<input HISTORYfile>` is the name of the DL_POLY HISTORY file which is missing its header
`<output HISTORYfile>` is the name of the new DL_POLY HISTORY file to write

## Switches

No switches exist for **addheader**.

