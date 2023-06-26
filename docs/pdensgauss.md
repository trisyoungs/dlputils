---
title: pdensgauss
brief: Apply Gaussian smooth to probability density
taxonomy:
  category: docs
  classification: "Spatial Probability Densities"
docroot: /dlputils/docs
---

**pdensgauss** applies a Gaussian smooth to an existing probability density, writing a new probability density file in the process. Probability density files calculated from [a]**pdens**,dlputils,manual,pdens[/a] can be noisy, depending on the number of frames / molecules used in the underlying averaging. Applying a Gaussian smooth can help to eliminate voxel noise, as well as smooth larger regions.

The applied function is a typical 3-dimensional Gaussian of the form:

<figure>
  <image>img/pdensgauss-gauss3.png</image>
</figure>

Related programs:
To calculate probability densities, see [**pdens**](/dlputils/docs/utilities/pdens).

## Basic Usage

```
pdensgauss <PDENSfile> <newPDENSfile> <sigma>
```

Where:
+ `<PDENSfile>` is the name of the existing PDENS file
+ `<newPDENSile>` is the name of the new PDENS file to be generated
+ `<sigma>` is a parameter controlling the width of the Gaussian used in the smoothing

## Switches

`-min` _value_

Prune the input pdens grid, removing all voxel values below the specified value.

`-max` _value_

Prune the input pdens grid, removing all voxel values above the specified value.

## Usage

The `sigma` parameter controls the width of the Gaussian used for smoothing the data in the probability density - a value of around 0.4 - 0.5 produces acceptable smoothing without removing too many subtleties in the features of the data. Note that the probability density is not treated as periodic when performing the smooth (i.e. there is assumed to be no data beyond the limits of the defined voxel grid).


