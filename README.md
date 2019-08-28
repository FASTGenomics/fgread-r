[![Build Status](https://travis-ci.org/FASTGenomics/fgread-r.svg?branch=master)](https://fastgenomics.github.io/fgread-r/docs/)

# Description

This package implements convenience functions for loading data sets in the
[FASTGenomics][fg] [analysis][fg_analysis] environment.  The functions from this package
will let you list and load data sets for which the analysis was defined.

[fg]: https://beta.fastgenomics.org/webclient/
[fg_analysis]: https://beta.fastgenomics.org/webclient/searchPage/analyses

## Usage

### Listing data sets

To list the data sets simply call the `fgread::list_datasets()` function

``` R
dsets_list <- fgread::list_datasets()
```

The `dsets_list` would then contain the information about the location, format, title,
etc. about of each data set.

```
[[1]]
id: 1
title: Loom data set
format: Loom
path: ../tests/data/readers/dataset_0001

[[2]]
id: 2
title: Seurat Object data set
format: Seurat Object
path: ../tests/data/readers/dataset_0002
```

Note, that `fgread::list_datasets()` does not load any of the data sets.  It's purpose
is to get a list of available data sets, from which you can select the ones you would
like to load.

### Loading a single data set

To load a single data set use `fgread::read_dataset`.  The code below loads the first
data set from the list (the "Loom data set") and returns a [Seurat][seurat] object

``` R
seurat <- fgread::read_dataset(dsets_list[[1]])
```

To load the second data set simply run

``` R
seurat2 <- fgread::read_dataset(dsets_list[[2]])
```

The `fgread::read_dataset` function resolves the underlying format of the data set
automatically, based on the `format` attributes contained in the `dsets_list[[1]]`.

[seurat]: https://satijalab.org/seurat/

### Loading multiple data sets

Similarly, one can load multiple data sets with a single command:
`fgread::read_datasets` (note the `s` at the end).  The command loads all available data
sets into _separate_ Seurat objects and returns a list of these objects (where the
indices correspond to the indices from `fgread::list_datasets`).

``` R
dsets <- fgread::read_datasets(dsets_list)
```
Now the `dsets` is a list containing two Seurat Objects

```
[[1]]
An object of class Seurat
16892 features across 298 samples within 1 assay
Active assay: RNA (16892 features)

[[2]]
An object of class Seurat
33538 features across 1222 samples within 1 assay
Active assay: RNA (33538 features)
```

Used without any arguments `fgread::read_datasets()` loads all data sets

``` R
dsets <- fgread::read_datasets()
```


```
[[1]]
An object of class Seurat
16892 features across 298 samples within 1 assay
Active assay: RNA (16892 features)

[[2]]
An object of class Seurat
33538 features across 1222 samples within 1 assay
Active assay: RNA (33538 features)
```

# Supported formats

The following formats are supported by this package
- [AnnData](https://github.com/theislab/anndata)
- [CellRanger (hdf5)](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices)
- [Drop-Seq (tsv)](https://github.com/Hoohm/dropSeqPipe/)
- [Loom](http://loompy.org/)
- [Seurat Object](https://satijalab.org/seurat/)

# Known issues

- Loading a relatively modest Drop-Seq data set (20k cell barcodes) uses around 10GB
  peak memory.  This could go over the 16GB limit with larger data sets.

- There is a function for reading AnnData in Seurat but it's buggy and does not seem to
  work on some test data sets.  Perhaps this will be fixed in future releases of Seurat
  but for now we use a custom implementation that only reads `.X`, `.obs` and `.var`
  components.  The function is also limited to count tables in the CSR format.

# Development and testing

Clone the repository along with the test data by running

``` bash
git clone git@github.com:FASTGenomics/fgread-r.git
cd fgread-r
git submodule init
git submodule update
```

Run R the `fgread-r` directory, install the devtools package (if you don't have it already)

``` R
install.packages("devtools")
```

And install the package dependencies

``` R
devtools::install_deps(upgrade="never")
```

Once the dependencies are there you can run the tests with

``` R
devtools::test()
```
