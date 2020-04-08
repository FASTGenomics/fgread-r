context("Listing datasets using get_datasets()")

DATADIR <- "../data/all_datasets/"

test_that("Check the list function", {
    dsets_list <- fgread::get_datasets(data_dir=DATADIR)
    dirs <- list.dirs(path = DATADIR, recursive = F)
    expect_equal(length(dsets_list), length(dirs))
})
