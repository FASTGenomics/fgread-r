context("Listing data sets")

DATADIR <- "../data/readers/"

test_that("Check the list function", {
    dsets_list <- fgread::list_datasets(data_dir=DATADIR)
    expect_equal(length(dsets_list), 5)
})
