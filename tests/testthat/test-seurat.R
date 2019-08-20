context("Loading into Seurat")

DATADIR <- "../data/readers"

## list
DATASETS <- list(
    list(id=1, dim=c(16892, 298), title="Loom data set", format="Loom"),
    list(id=2, dim=c(33538, 1222), title="Seurat Object data set", format="Seurat Object"),
    list(id=3, dim=c(20, 10), title="AnnData data set", format="AnnData"),
    list(id=4, dim=c(33538, 1222), title="10x (hdf5) data set", format="10x (hdf5)"),
    list(id=5, dim=c(99, 20000), title="Drop-Seq (tsv) data set", format="Drop-Seq (tsv)")
)

Map(
    function(dset){
        test_that(glue::glue("Format {dset$format} loads"),
                  suppressWarnings({
                      dsets_list <- fgread::list_datasets(data_dir=DATADIR)
                      seurat <- fgread::read_dataset(dsets_list[[dset$id]])
                      expect_equal(dim(seurat), dset$dim)
                      expect_equal(seurat@project.name, dset$title)
                      expect_equal(seurat@misc$fg_metadata, dsets_list[[dset$id]]@metadata)
                  }))
    },
    DATASETS
)
