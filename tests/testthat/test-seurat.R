context("Loading into Seurat using load_data()")

DATADIR <- "../data/all_datasets"
readers <- DEFAULT_READERS
supported_readers_str <- paste(names(readers), collapse = ", ") #TODO: Test multi expression files explicit

test_that("Check equality of metadata in json and Seurat object", {
  for (dir in list.dirs(path = DATADIR, recursive = F)) {
    json_path <- file.path(dir, INFO_FILE_NAME)
    json_info <- jsonlite::read_json(json_path)
    print(json_info$title)
    file <- json_info["expressionDataFileInfos"][[1]][[1]]$name
    format <- tail(strsplit(file, "\\.")[[1]], n = 1)
    if (format %in% names(readers)) {
      # run tests with and without experimental readers
      seurat <- tryCatch({
        fgread::load_data(json_info$title, experimental_readers = T, data_dir = DATADIR)
      }, error = function(e) {
        fgread::load_data(json_info$title, experimental_readers = F, data_dir = DATADIR)
      })
      expect_equal(dim(seurat)[[1]], json_info$numberOfGenes)
      expect_equal(dim(seurat)[[2]], json_info$numberOfCells)
      expect_equal(seurat@project.name, json_info$title)
      expect_equal(seurat@misc$fastgenomics$id, json_info$id)
    }
  }
})
