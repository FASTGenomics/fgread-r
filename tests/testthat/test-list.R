context("Listing datasets")

DATADIR <- "../data/all_datasets/"

test_that("Check the list function", {
    dsets_list <- fgread::get_datasets(data_dir=DATADIR)
    expect_equal(length(dsets_list), 11)
})
