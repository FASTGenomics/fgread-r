[![Build Status](https://travis-ci.org/FASTGenomics/fgread-r.svg?branch=master)](https://fastgenomics.github.io/fgread-r/docs/)

# Description

This package implements convenience functions for loading datasets in the
[FASTGenomics][fg] [analysis][fg_analysis] environment.  The functions from this package
will let you list and load datasets for which the analysis was defined.

[fg]: https://beta.fastgenomics.org/webclient/
[fg_analysis]: https://beta.fastgenomics.org/webclient/searchPage/analyses

## Usage

### Listing datasets

To list the datasets simply call the `fgread::get_datasets()` function

``` R
dsets_list <- fgread::get_datasets()
```

The `dsets_list` would then contain the information about the location, format, title,
etc. about of each dataset.

```
[[1]]
id: 1
title: Loom dataset
format: Loom
path: ../tests/data/readers/dataset_0001

[[2]]
id: 2
title: Seurat Object dataset
format: Seurat Object
path: ../tests/data/readers/dataset_0002
```

Note, that `fgread::get_datasets()` does not load any of the datasets.  It's purpose
is to get a list of available datasets, from which you can select the ones you would
like to load.

### Loading a single dataset

To load a single dataset use `fgread::read_dataset`.  The code below loads the first
dataset from the list (the "Loom dataset") and returns a [Seurat][seurat] object

``` R
seurat <- fgread::read_dataset(dsets_list[[1]])
```

To load the second dataset simply run

``` R
seurat2 <- fgread::read_dataset(dsets_list[[2]])
```

The `fgread::read_dataset` function resolves the underlying format of the dataset
automatically, based on the `format` attributes contained in the `dsets_list[[1]]`.

[seurat]: https://satijalab.org/seurat/

### Loading multiple datasets

Similarly, one can load multiple datasets with a single command:
`fgread::read_datasets` (note the `s` at the end).  The command loads all available data
sets into _separate_ Seurat objects and returns a list of these objects (where the
indices correspond to the indices from `fgread::get_datasets`).

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

Used without any arguments `fgread::read_datasets()` loads all datasets

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

- Loading a relatively modest Drop-Seq dataset (20k cell barcodes) uses around 10GB
  peak memory.  This could go over the 16GB limit with larger datasets.

- There is a function for reading AnnData in Seurat but it's buggy and does not seem to
  work on some test datasets.  Perhaps this will be fixed in future releases of Seurat
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
