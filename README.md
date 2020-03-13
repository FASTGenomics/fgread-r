[![Build Status](https://travis-ci.org/FASTGenomics/fgread-r.svg?branch=master)](https://fastgenomics.github.io/fgread-r/docs/)

# FASTGenomics Reader Module for R

This package implements convenience functions for loading datasets in the
[FASTGenomics][fg] analysis environments. The functions from this package
let you list and load datasets associated to the corresponding analysis.

[fg]: https://beta.fastgenomics.org/

## Documentation

For the general documentation on how to use the reader, please visit our FASTGenomics Documentation.

For details on the available functions see the [API Documentation](https://fastgenomics.github.io/fgread-r/docs/).

## Known issues

- Seurat's loading function uses a lot of memory (already 10GB RAM for a relatively modest Drop-Seq dataset 
  of 20k cell barcodes). For larger datasets one may therefore exceed the resource limits on FASTGenomics.

- There is no functionality to read AnnData objects in Seurat. Therefore, we provide a custom reader. 
  This custom reader, however, only reads the `.X`, `.obs` and `.var` components of the AnnData object.
  The function is also limited to count tables in the CSR format.

## Development and testing

Clone the repository along with the test data by running

```bash
git clone git@github.com:FASTGenomics/fgread-r.git
cd fgread-r
git submodule init
git submodule update
```

Run R in the `fgread-r` directory, install the devtools package (if you don't have it already)

```R
install.packages("devtools")
```

and install the package dependencies

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
