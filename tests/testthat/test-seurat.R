context("Loading into Seurat")

DATADIR <- "../data/readers"

## list
DATASETS <- list(
    list(id=1, dim=c(16892, 298), title="Loom dataset", format="Loom"),
    list(id=2, dim=c(33538, 1222), title="Seurat Object dataset", format="Seurat Object"),
    list(id=3, dim=c(20, 10), title="AnnData dataset", format="AnnData"),
    list(id=4, dim=c(33538, 1222), title="10x (hdf5) dataset", format="10x (hdf5)"),
    list(id=5, dim=c(99, 20000), title="tab-separated text dataset", format="tab-separated text"),
    list(id=10, dim=c(99, 499), title="comma-separated text dataset", format="comma-separated text"),
    list(id=11, dim=c(99, 499), title="tab-separated text variant dataset", format="tab-separated text")
)

Map(
    function(dset){
        test_that(glue::glue("Format {dset$format} loads"),
                  suppressWarnings({
                      dsets_list <- fgread::get_datasets(data_dir=DATADIR)
                      seurat <- fgread::read_dataset(dsets_list[[dset$id]])
                      expect_equal(dim(seurat), dset$dim)
                      expect_equal(seurat@project.name, dset$title)
                      expect_equal(seurat@misc$fastgenomics$metadata, dsets_list[[dset$id]]@metadata)
                      expect_equal(seurat@misc$fastgenomics$id, dset$id)
                  }))
    },
    DATASETS
)
