[![Build Status](https://travis-ci.org/FASTGenomics/fgread-r.svg?branch=master)](https://fastgenomics.github.io/fgread-r/docs/)

# FASTGenomics Reader Module for R

This package implements convenience functions for loading datasets in the
[FASTGenomics][fg] [analysis][fg_analysis] environment. The functions from this package
will let you list and load datasets for which the analysis was defined.

[fg]: https://beta.fastgenomics.org/webclient/
[fg_analysis]: https://beta.fastgenomics.org/webclient/searchPage/analyses

## Documentation

For the general documentation on how to use the reader, please visit our FASTGenomics Documentation.

For details on the available functions see the [API Documentation](https://fastgenomics.github.io/fgread-r/docs/).

## Known issues

- Loading a relatively modest Drop-Seq dataset (20k cell barcodes) uses around 10GB
  peak memory. This could go over the 16GB limit with larger datasets.

- There is a function for reading AnnData in Seurat but it's buggy and does not seem to
  work on some test datasets. Perhaps this will be fixed in future releases of Seurat
  but for now we use a custom implementation that only reads `.X`, `.obs` and `.var`
  components. The function is also limited to count tables in the CSR format.

## Development and testing

Clone the repository along with the test data by running

```bash
git clone git@github.com:FASTGenomics/fgread-r.git
cd fgread-r
git submodule init
git submodule update
```

Run R the `fgread-r` directory, install the devtools package (if you don't have it already)

```R
install.packages("devtools")
```

And install the package dependencies

```R
devtools::install_deps(upgrade="never")
```

Once the dependencies are there you can run the tests with

```R
devtools::test()
```

### Build Documentaion

```R
devtools::document()
```
