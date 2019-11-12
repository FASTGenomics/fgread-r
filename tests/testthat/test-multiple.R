context("Loading multiple datasets")

DATADIR <- "../data/all_datasets"

test_that("Loading multiple datasets", suppressWarnings({
    dsets_list <- fgread::get_datasets(data_dir=DATADIR)
    sel <- c(1,2,3)
    dsets <- fgread::read_datasets(dsets_list[sel], experimental_readers=T)

    expect_equal(length(dsets), length(sel))
    print(names(dsets))
    for (id in sel){
        expect_equal(dsets[[id]]@misc$fastgenomics$id, id)
    }
}))
