---
title: rdf
brief: Calculate centre-of-mass radial distribution functions
taxonomy:
  category: docs
docroot: /dlputils/docs
template: manpage
---

**rdf** calculates, by default, the centre-of-mass radial distribution functions between all species in the system, although the options allow a point other than the centre-of-mass to be used in the calculation. In addition, for surrounding molecules a different point can be defined to the central molecule.

Related programs:
+ To calculate individual atomic partials for one species, see [**rdfaa**](/dlputils/docs/rdfaa).
+ To calculate site-COM partials, see [**rdfss**](/dlputils/docs/rdfss).
+ To calculate site-site partials dependent on other sites, see [**rdfdep**](/dlputils/docs/rdfdep).

## Basic Usage

```
rdf <HISTORYfile> <OUTPUTfile>
```

Where:
`<HISTORYfile>` is the name of the DL_POLY HISTORY file
`<OUTPUTfile>` is the name of the DL_POLY OUTPUT file

## Switches

`-bin` _width_
Specify the histogram bin _width_ to use in the calculation (default = 0.1).

`-compair` _sp_ _i_ _j_
Instead of the centre-of-mass, use the average coordinates of atom indices _i_ and _j_ for species index _sp_. The atom indices are local atom indices for the species in question (i.e. range from 1 to the number of atoms in the molecule).

`-discard` _n_
Discard _n_ frames at the start of the trajectory.

`-frames` _n_
Stop calculating after processing _n_ frames.

`-header` _file_
Read in trajectory header from the specified DL_POLY HIStory _file_. Useful when the desired target frames are present in a restarted (and not appended) HIStory file.

`-nonorm`
Do not normalise the radial distribution functions to the bulk number density of species.

`-otherpair` _sp_ _i_ _j_
Instead of the centre-of-mass (or `-compair` specification), use the average coordinates of atom indices _i_ and _j_ for species index _sp_ when considering the position of surrounding molecules. The atom indices are local atom indices for the species in question (i.e. range from 1 to the number of atoms in the molecule).

`-zminus` _z_
Surrounding molecules will only be accepted if the z-value of it's centre-of-mass (or defined point) is less than that of the central molecule. Mutually exclusive with `-zplus`.

`-zplus` _z_
Surrounding molecules will only be accepted if the z-value of it's centre-of-mass (or defined point) is greater than that of the central molecule. Mutually exclusive with `-zminus`.

## Normalisation

All RDFs are normalised to reflect the number density of the surrounding species, as is standard practice. However, for RDFs where the central molecule species is the same as the surrounding one, the number density is calculated using N-1 molecules.

## Output

Given an input HISTORY file `results.HISu`, then output files are as follows:
`results.rdfNNMM` : RDF between species NN and MM


