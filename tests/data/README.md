# Test datasets for FASTGenomics

This folder contains small datasets with all currently supported formats. It is intended to be used for testing purposes and during development. A testing environment is provided in the git repository [analysis-test-environment](https://github.com/FASTGenomics/analysis-test-environment).

## Structure

* datasets: folder with symlinks to all_datasets with naming according to data format
  * 10x_h5: 10x (hdf5) format, synthetic test data [1222 cells / 33538 genes]
  * 10x_mtx: 10x mtx legacy, synthetic test data [30 cells / 1000 genes]
  * 10x_mtx_v3: 10x mtx v3, synthetic test data [30 cells / 1000 genes]
  * anndata_h5ad: Anndata, synthetic test data [10 cells / 20 genes]
  * csv: dense matrix, synthetic test data [499 cells / 99 genes]
  * loom: synthetic test data [298 cells / 16892 genes]
  * notset_no_exist: data format "not set", for testing purpose only. Data file missing
  * other_no_exist: data format "other", for testing purpose only. Data file missing
  * seurat_rds: Seurat object, synthetic test data [1222 cells / 33538 genes]
  * tsv: dense matrix, first element in header missing ("gene"), synthetic test data [20000 cells / 99 genes]
  * tsv_var: dense matrix, synthetic test data [499 cells / 99 genes]
* all_datasets: Small test datasets with all currently supported formats in FASTGenomics named according to FASTGenomics terminology, i.e., datasets_xxxx
* README.md: this readme file
