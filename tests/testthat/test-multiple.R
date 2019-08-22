context("Loading multiple data sets")

DATADIR <- "../data/readers"

test_that("Loading multiple data sets", suppressWarnings({
    dsets_list <- fgread::list_datasets(data_dir=DATADIR)
    sel <- c(1,2,3)
    dsets <- fgread::read_datasets(dsets_list[sel])

    expect_equal(length(dsets), length(sel))
    print(names(dsets))
    for (id in sel){
        expect_equal(dsets[[id]]@misc$fastgenomics$id, id)
    }
}))
