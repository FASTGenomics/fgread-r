context("Loading into Seurat using load_data()")

DATADIR <- "../data/all_datasets"

## list
DATASETS <- list(
    list(id = "dataset-f93640431b0e4835a55332mbrnepxlne", dim = c(16892, 298), title = "Loom dataset", format = "Loom"),
    list(id = "dataset-f93640431b0e4835a55332r9ii3nzffn", dim = c(33538, 1222), title = "Seurat Object dataset", format = "Seurat Object"),
    list(id = "dataset-f93640431b0e4835a55332figaqdiem3", dim = c(20, 10), title = "AnnData dataset", format = "AnnData"),
    list(id = "dataset-f93640431b0e4835a55332bruktcoddp", dim = c(33538, 1222), title = "10x (hdf5) dataset", format = "10x (hdf5)"),
    list(id = "dataset-f93640431b0e4835a5533200fet8a6jm", dim = c(99, 20000), title = "tab-separated text dataset", format = "tab-separated text"),
    list(id = "dataset-f93640431b0e4835a55332ssdhigwux7", dim = c(1000, 30), title = "mtx legacy dataset", format = "10x (mtx)"),
    list(id = "dataset-f93640431b0e4835a553327j4gmcofeu", dim = c(1000, 30), title = "mtx v3 dataset", format = "10x (mtx)"),
    list(id = "dataset-f93640431b0e4835a55332k8vohrjpyv", dim = c(99, 499), title = "comma-separated text dataset", format = "comma-separated text"),
    list(id = "dataset-f93640431b0e4835a55332pj5jogcig7", dim = c(99, 499), title = "tab-separated text variant dataset", format = "tab-separated text"),
    list(id = "dataset-f93640431b0e4835a55332fpgkljaqdz", dim = c(20, 10), title = "AnnData dense dataset", format = "AnnData")
)


Map(
    function(dset) {
      test_that(glue::glue("Format {dset$format} loads"),
                  suppressWarnings({
        data_table <- fgread::ds_info(data_dir = DATADIR, pretty = FALSE)
        print(dset$title)
        seurat <- fgread::load_data(dset$title, experimental_readers = T, data_dir = DATADIR)
        expect_equal(dim(seurat), dset$dim)
        expect_equal(seurat@project.name, dset$title)
        # expect_equal(seurat@misc$fastgenomics$metadata, data_table[[dset$id]]@metadata)
        expect_equal(seurat@misc$fastgenomics$id, dset$id)
      }))
    },
    DATASETS
)
