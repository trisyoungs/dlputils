---
title: Usage
brief: General usage notes
taxonomy:
  category: docs
docroot: /dlputils/docs
header_class: alt
---

Most of the codes in **dlputils** package have a similar usage syntax on the command line - at the simplest level, the names of DL_POLY HIStory and OUTput files are the only required arguments. For example:

```
user@pc:~> rdf mysystem.HISu mysystem.OUT
```

Note that the extension of the HIStory file *is* important - passing a default-named HISTORY file makes it difficult to detect whether the file is formatted or unformatted. For this reason, all the **dlputils** assume that a trajectory with the extension _.HISu_ is unformatted, while one with _.HISf_ is formatted. The extension of the OUTput file is not important. It is not necessary to specify the number of molecule types (species) defined in the system, number of atoms in molecules, etc., since this is read from the supplied OUTput file. You can also supply a _0_ (zero) in place of the OUTput file name, in which case you will be asked for the number of molecular species, and the correct number of molecules / atoms for each. The codes are designed to read DL_POLY v2 and DL_CLASSIC HIStory and OUTput files, but making them compatible for DL_POLY v3 (and beyond) should not be a difficult task (modify dlprw as necessary).

Since the first two arguments are mandatory, running any of the codes with no arguments will print out usage information, including a list of available command-line switches which control various relevant parameters in the codes.

In most cases one of more output files are generated, rather than results being written to screen.  The names of these files is always based on the name of the HIStory file, but only if the filename has an extension (it doesn't matter what the extension is). This is removed and replaced by a new, descriptive extension. For example, the **rdf** program produces files called _.rdfNM_, where _N_ and _M_ are the identifying species numbers. So, the above example will produce output files called _mysystem.rdf11_, _mysystem.rdf12_, etc., depending on the number of species in the system.

##Not Using DL_POLY?

Although specific read-write routines for DL_POLY trjectories are the only ones currently implemented, that does not necessarily mean you can't use **dlputils** to analyse your data. As long as you know the unit cell parameters and a standard xyz trajectory, this can be converted to a DL_POLY history file. A suitable OUT file (which is used to describe the frame contents in terms of species, molecules, and atoms) can also be constructed quite easily. See the [a]**xyz2his**,dlputils,manual,xyz2his[/a] [a]manual page,dlputils,manual,xyz2his[/a] for specific instructions.

