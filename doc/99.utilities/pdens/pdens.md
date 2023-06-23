---
title: pdens
brief: Calculate three-dimensional (spatial) probability densities
visible: true
taxonomy:
  category: docs
  classification: "Spatial Probability Densities"
docroot: /dlputils/docs
template: manpage
---

**pdens** calculates the spatial probability densities of molecules around other molecules, and amounts to performing 3-dimensional histogram binning (into volume elements or 'voxels') as opposed to 1-dimensional histogram binning in the case of the radial distribution function (averaging the spatial probability density spherically gives the standard radial distribution function). **pdens** uses the defined axis origins to represent the positions of the surrounding molecules or, if one has not beeb given, the centre-of-mass is used instead. Different points on the surrounding molecules may be defined with the `-otherpair` option.


For the central species, a 3D histogram is constructed for each molecular species in the system, extending a finite number of voxels of a given size in each of the Cartesian axes (controlled by the `-delta` and `-grid` keywords respectively). The intramolecular spatial probabilty density is also calculated for the central species, and provides information on the conformational space explored by the central molecules (and 'constrained' by the defined axis system).


If the probability of finding specific orientations of surrounding molecules is needed, this can be achieved by using the `-orient` keyword. This permits surrounding molecules to be excluded based on the angular difference between a specific axis on the central and surrounding molecules.


Note that, since binning is performed in three dimensions, many more frames are required compared to [**rdf**](/dlputils/docs/utilities/rdf) in order to get good statistics (and hence better-looking surfaces). A post-processing utility to perform Gaussian smoothing of the data, [**pdensgauss**](/dlputils/docs/utilities/pdensgauss), is available.

Related programs:
+ To apply a Gaussian smooth to an existing probability density, see [**pdensgauss**](/dlputils/docs/utilities/pdensgauss).

## Basic Usage

```
pdens <HISTORYfile> <OUTPUTfile> <sp> <othersp>
```

Where:
+ + `<HISTORYfile>` is the name of the DL_POLY HISTORY file
+ + `<OUTPUTfile>` is the name of the DL_POLY OUTPUT file
+ + `<sp>` is the integer index of the species to consider as the central one, around which the probability densities will be calculated
+ + `<othersp>` is a comma-separated list of indices specifying which species to consider as the surrounding ones

Although this is enough to run **pdens**, axis definitions for at least the central species must be provided (see the section on [**Defining Axes**](/dlputils/docs/utilities/pdens#axes) below).

## Switches

`-atoms` _sp_ _list_


For species _sp_, use the coordinates of each individual atom index in the supplied (and comma-separated) _list_ as a point in the surrounding species, with each contributing individually to the final probability density. These coordinates will be used in preference to the coordinate centre derived from -axis and -otheraxis specifications. This option is mutually-exclusive with -cog.

`-axis` _sp_ _x1_ _x2_ _y1_ _y2_


Define the axis system for species index _sp_ from the four supplied atom indices. The atom indices are local atom indices for the species in question (i.e. range from 1 to the number of atoms in the molecule). See the section on [**Defining Axes**](/dlputils/docs/utilities/pdens#axes) for more information.

`-cartesian` _sp_


For the specified species, use basic Cartesian axes (i.e. a set of axes invariant with molecular orientation).

`-cog` _sp_ _list_


For species _sp_, use the average centre-of-geometry of the atom indices in the supplied (and comma-separated) _list_ as the point to consider in the surrounding species. This centre will be used in preference to the coordinate centre derived from -axis and -otheraxis specifications. This option is mutually-exclusive with -atoms.

`-delta` _d_


Set the size of the intermolecular histogram bin voxel to be _d_ along each side (default = 0.5).

`-end` _n_


Stop calculating when frame _n_ is reached.

`-grid` _npoints_


Specify that the intermolecular voxel histogram should extend _npoints_ in each positive and negative Cartesian direction (default = 20). Note that the overall extent of the volume covered also depends on the voxel size (see the `-delta` option).


`-header` _file_


Read in trajectory header from the specified DL_POLY HIStory _file_. Useful when the desired target frames are present in a restarted (and not appended) HIStory file.

`-intra`


Turn on calculation of the intramolecular probability distribution.

`-maxdist` _r_


Set tha maximum allowable distance between the central molecule and the surrounding molecule to _r_.

`-mindist` _r_


Set tha minimum allowable distance between the central molecule and the surrounding molecule to _r_.

`-molmap` _file_


Instead of considering all molecules of the central species from each frame, the supplied _file_ dictates specific indices of the central molecule in each frame that should be considered.

`-nointer`


Turn off calculation of the intermolecular probability distributions.

`-orient` _sp_ _axis_ _angle_ _delta_


By default all surrounding molecules are binned, regardless of their orientation. The `-orient` keyword allows deviations in angle between the central and surrounding molecule to be defined for each axis, outside of which the molecules are omitted. See the section on [**Selecting Molecule Orientations**](/dlputils/docs/utilities/pdens#orientations) for more information.


`-otheraxis` _sp_ _x1_ _x2_ _y1_ _y2_


The `-otheraxis` keyword allows a separate set of axes to be defined for a given species which will be used to define the centre and orientation of surrounding molecules, distinct from the axes defining the central molecule. These axes, if specified, will be used in orientation calculations rather than those provided by -axis. Note that if either -atoms or -cog is given for the same species, those will be used in preference to the axes centre inferred here.


`-pdelta` _d_


Set the size of the intramolecular histogram bin voxel to be _d_ along each side (default = 0.15).

`-pgrid` _npoints_


Specify that the intramolecular voxel histogram should extend _npoints_ in each positive and negative Cartesian direction (default = 50). Note that the overall extent of the volume covered also depends on the voxel size (see the `-pdelta` option).


`-select` _sp_ _pdensfile_ _atoms_ _cutoff_


Restricts how the current pdens for species _sp_ is constructed, by making it dependent on the pdens of a previous calculation on the same species. For the target species _sp_, the existing _pdensfile_ is loaded and used as a reference. The point used to calculate this existing density must be specified using the comma-separated list of atom indices in _atoms_ - the geometric mean of these atom coordinates is used. When calculating the new probability density, a point is added to the average only if the corresponding point specified by _atoms_ in the existing density is equal to or above the _cutoff_ specified. In this manner, a distribution of another point on a molecule, for instance, can be constructed relative to the probabilty of finding another point on the same molecule at some density.

`-start` _n_


Start calculating on trajectory frame _n_ (default = 1).

## Output <a id="output"></a>

Given an input HISTORY file **results.HISu**, then output files are as follows:
+ **results.NNMM.pdens** : probability density of point on species MM around the central species NN
+ **results.intraNN.pdens** : intramolecular probability density of species NN
+ **results.avgNN.xyz** : average coordinates of the species NN

## Axis Definitions <a id="axes"></a>

For any calculation of the spatial density around a central molecule, an axes definition must be provided that uniquely defines the local orientation of the molecule (or at least part of it). This means that the central species must contain at least three atoms in a non-linear arrangement. The exception to this rule is when using Cartesian (invariant) axes for the central molecule (see the `-cartesian` keyword above), which will work even for atomic species.


Four atom indices are provided in order to define a set of axes - the first two, _x1_ and _x2_ define the direction of the x axis and also the origin of the axis system, taken to be the average of the coordinates of the two atoms. Because of this, the two indices _x1_ and _x2_ [b]must[/b] be different. The y axis is then formed from the axis origin to the average coordinates of the two atom indices _y1_ and _y2_. Note that _y1_ and _y2_ may be the same atom index, and also that they do not have to represent [i]exactly[/i] the direction of the y axis since the vector is orthogonalised with respect to the already defined x axis. The z axis is then formed from the cross product of the x and y axes.

The atoms used to define the axis system effectively 'fixes' those atoms as a point of reference for the central molecule - unless the molecule is quite rigid, the other atoms may 'average out' into a tangled ball in the corresponding intramolecular spatial density functions. Furthermore, for non-rigid molecules, the further one moves away from the atoms used to define the axis system, the more 'diffuse' and averaged the spatial density is likely to become. As such, for molecules that possess more than one site or functional group of interest it is prudent to run separate calculations with different axis definitions. For instance, for the methanol molecule we may wish to focus on both the methyl and hydroxyl ends of the molecule, and so might make the following axis definitions (performing a separate **pdens** calculation for each):

![](pdens-axes.png){.imgfull}

This shows two possible axis definitions for a methanol molecule, with the middle image focussing on the hydroxyl group (`-axis 1 1 3 2 2`) and the right-hand image focussing on the methyl group (`-axis 1 6 5 4 4`).


## Selecting Molecule Orientations <a id="orientations"></a>

By default all molecules contribute to the spatial density, provided they are within the range of the defined grid. In some circumstances, however, it is desirable to construct the spatial density of only those molecules with specific orientations relative to the central one. In **pdens** this can be achieved by defining one or more angle deviations between the axes of the central and surrounding molecule, outside of which the second molecule will be ignored.

Consider the case of benzene, where we have defined the axes as follows:

![](pdens-benzene-axes.png){.imgfull}

As an example, let us suppose that we wish to find out in which positions co-planar molecules sit around a central molecule - in other words, how do molecules sit parallel to each other in close proximity. Since the z-axis is perpendicular to the plane of the benzene ring, we can specify an angular deviation away from this axis to request only parallel molecules are considered. Within **pdens**, this amounts to considering the dot product of the z axis of the central molecule with that of the candidate molecule:

![](pdens-benzene-orient.png){.imgfull}

We may request that only such parallel orientations are considered by passing `-orient 1 3 0.0 10.0`, which would ask that species 1 (benzene) molecules should only be considered if the angle between the third axes (the z axes) is within 10.0 degrees of 0.0. However, the present example is actually a special case, since the XY plane is also a plane of symmetry, and so when the angle theta tends towards 180 degrees this also corresponds to a parallel orientation. In such circumstances a negative _delta_ can be given, which basically indicates that a mirror plane exists for the axis in question.



