---
title: xyz2his
brief: Convert a sequence of xyz trajectory frames to an unformatted DL_POLY HISTORY file
visible: true
taxonomy:
  category: docs
  classification: "History File Conversion / Manipulation"
docroot: /dlputils/docs
template: manpage
---

The use of **dlputils** to analyse output from other codes is possible, provided that the relevant trajectory information is converted to DL_POLY format first. Probably the easiest way is to get your data of interest written as an appended sequence of xyz format frames in a single file, and then use **xyz2his** to convert the data.

## Basic Usage

```
xyz2his <XYZfile> <HISTORYfile>
```

Where:
+ `<XYZfile>` is the name of the XYZ file containing the source frames
+ `<HISTORYfile>` is the target name of the DL_POLY HISTORY file to be written

## Switches

`-cell` _ax_ _ay_ _az_ _bx_ _by_ _bz_ _cx_ _cy_ _cz_

Specify the unit cell of the system by providing a full 3x3 matrix. Either this or `-cubic` should be used.


`-cubic` _a_

Specify a cubic unit cell for the system by providing the side length _a_. Either this or `-cell` should be used.


`-extraline`

Indicates that the xyz trajectory has a blank line between consecutive frames.

## XYZ File Format
Firstly, the expected format of the xyz file is as follows, using the simple case of a single isolated water molecule as an example:

```
3
Water - Frame 1
O            -0.000000     0.000000    -0.000000
H             0.742176     0.021779     0.590470
H            -0.786228     0.006122     0.530414
3
Water - Frame 2
O            -0.141170    -0.094786    -0.005418
H             0.746772     0.022545     0.590150
...
```

The first line of each frame contains the number of atoms to follow, and which must be the same for each and every subsequent frame. The second line contains the title of the frame, but may simply be a blank line if necessary. What follow is then the atom elements/names and their coordinates, one per line. The spacing between the data on the lines is not important. This is the standard format expected by the **xyz2his** program, but some programs write an extra blank line between frames (EPSR for example). If this is the case, you must give the `-extraline` option to **xyz2his**.


## Unit Cell Specification

It is assumed that the primary purpose of **dlputils** is for the analysis of data from simulations of periodic systems - as such a unit cell must be provided (since the DL_POLY HISTORY file format contains one). Either the full 3x3 unit cell matrix can be specified with the `-cell` option, or a cubic cell can be specified with the `-cubic` option and the unit cell length.


## Performing the Conversion

Assuming the input trajectory is _trajectory.xyz_, then run **xyz2his** as follows:

_With a full 3x3 cell specification..._
```
user@pc:~> xyz2his trajectory.xyz trajectory.HISu -cell 14.0 0.0 0.0 -4.0 34.0 0.0 -0.2 -1.8 10.0
```

_With a cubic cell specification..._
```
user@pc:~> xyz2his trajectory.xyz trajectory.HISu -cubic 37.5
```

## Creating an Accompanying OUT File

As stated earlier, most every program in **dlputils** requires that the DL_POLY OUT file be given as well as the HISTORY file, in order that the number of molecular species, as well as the number of molecules of each and the number of constituent atoms, may be determined automatically. In fact, only a few key lines are searched for in the OUT file, and so it is perfectly possible to construct one when one does no exist (when DL_POLY was not used, for instance). The most basic OUT file for a simple simulation of water is as follows:

```
SYSTEM SPECIFICATION

 number of molecular types               1

 molecular species type                  1
 name of species:             Water
 number of molecules                  300
 number of atoms/sites                  3
```

It is important to note that, while the number of spaces between the descriptive text at the start of each line and the property that follows is unimportant, the single spaces at the beginning of each line (except the one beginning 'SYSTEM SPECIFICATION') are expected and necessary, since this is how DL_POLY writes it. The content is self-explanatory - the number of different molecule types is specified, and then there are three lines for each molecule type, giving the name, number of molecules, and number of atoms per molecule. This information is repeated for each molecule type in the system, so for a water/methanol mixture the required OUT file becomes:

```
SYSTEM SPECIFICATION

number of molecular types               2

molecular species type                  1
name of species:             Water
number of molecules                  300
number of atoms/sites                  3

molecular species type                  2
name of species:             Methanol
number of molecules                  150
number of atoms/sites                  6
```
