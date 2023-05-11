test_that("block preprocessing by number of columns", {
    prot_num_cols <- prot / ncol(prot)

    expect_equal(block_preproc(prot, "num_cols"), prot_num_cols)
})
