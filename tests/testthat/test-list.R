context("Listing datasets")

DATADIR <- "../data/readers/"

test_that("Check the list function", {
    dsets_list <- fgread::get_datasets(data_dir=DATADIR)
    expect_equal(length(dsets_list), 8)
})
