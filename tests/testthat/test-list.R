context("Listing datasets using ds_info()")

DATADIR <- "../data/all_datasets/"

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
