context("Listing datasets using ds_info()")

DATADIR <- "../data/all_datasets/"
INFO_FILE_NAME <- "dataset_info.json"

test_that("Check the info function", {
  data_table <- fgread::ds_info(data_dir = DATADIR)
  dirs <- list.dirs(path = DATADIR, recursive = F)
  expect_equal(dim(data_table)[1], length(dirs))
})

test_that("Check the info function for one dataset", {
  data_table <- fgread::ds_info('Loom dataset', data_dir = DATADIR, pretty = FALSE, output = TRUE)
  expect_equal(dim(data_table)[1], 1)
})

test_that("Check the info function for one dataset", {
  data_table <- fgread::ds_info('dataset-f93640431b0e4835a55332r9ii3nzffn', data_dir = DATADIR, pretty = FALSE, output = TRUE)
  expect_equal(dim(data_table)[1], 1)
})

test_that("Check equality of metadata in json and list", {
  for (dir in list.dirs(path = DATADIR, recursive = F)) {
    json_path <- file.path(dir, INFO_FILE_NAME)
    json_info <- jsonlite::read_json(json_path)
    json_info["schemaVersion"] <- NULL
    ds_info <- fgread::ds_info(json_info$title, output = T, pretty = F)
    for (col in names(json_info)) {
        expect_equivalent(json_info[col], ds_info[col][[1]])
    }
  }
})
