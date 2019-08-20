[![Build Status](https://travis-ci.org/FASTGenomics/fgread-r.svg?branch=master)](https://fastgenomics.github.io/fgread-r/docs/)

# Supported formats

- CellRanger (hdf5), via Seurat
- Loom, custom implementation
- AnnData, custom implementation
- Drop-Seq (tsv), custom implementation
- Seurat Object, via Seurat

# Known issues

- Loading a relatively modest Drop-Seq data set (20k cell barcodes) uses around 10GB
  peak memory.  This could go over the 16GB limit with larger data sets.

- There is a function for reading AnnData in Seurat but it's buggy and does not seem to
  work on my test data sets.  Perhaps this will be fixed in future releases of Seurat
  but for now we use a custom implementation that only reads `.X`, `.obs` and `.var`
  components.

- The current AnnData reader works only when the `.X` component is in the sparse CSR
  format.  Other formats would had to be implemented in the future.


# Development and testing

Clone the repository along with the test data by running

``` bash
git clone git@github.com:FASTGenomics/fgread-r.git
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
testthat::test()
```
