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

- For now, the file types are taken from a json file.  This will have to be updated once
  the proper architecture is in place.
