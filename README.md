[![Build Status](https://travis-ci.org/FASTGenomics/fgread-r.svg?branch=master)](https://fastgenomics.github.io/fgread-r/docs/)

# FASTGenomics Reader Module for R

This package implements convenience functions for loading datasets in the
[FASTGenomics][fg] [analysis][fg_analysis] environment. The functions from this package
will let you list and load datasets for which the analysis was defined.

[fg]: https://beta.fastgenomics.org
[fg_analysis]: https://beta.fastgenomics.org/webclient/searchPage/analyses

## Documentation

For the general documentation on how to use the reader, please visit our [FASTGenomics Documentation][docs].

For details on the available functions see the [API Documentation](https://fastgenomics.github.io/fgread-r/docs/).

[docs]: https://beta.fastgenomics.org/docs

## Limitations

### TSV/CSV Reader

Loading large files, especially if they are very sparse, in plain text format need a lot of memory and can be rather slow.
For large datasets we recommend to use formats that support sparse data.

### AnnData Reader

As the internal Seurat reader for AnnData fails to convert certain datasets (depending on scaling, normalization etc.), we provide a custom reader.
This custom reader, however, only reads the `.X`, `.obs` and `.var` components of the AnnData object.

## Known issues

Please report the issues through [github][issues].

[issues]: https://github.com/FASTGenomics/fgread-r/issues

## Development and testing

Clone the repository along with the test data by running

```bash
git clone --recurse-submodules git@github.com:FASTGenomics/fgread-r.git
```

Then enter the `fgread-r` directory and install the dependencies with

```R
install.packages("devtools")
devtools::install_deps(upgrade="never")
```

To test the package use

```R
devtools::test()
```

### Build Documentaion

```R
devtools::document()
```
