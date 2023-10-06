test_that("only existing omics", {
    expect_error(ord_loadings(mcia_out = mcia_results, omic = "test"))
})

# all omics and factor 1, ranked in descending order
all_pos_1 <- ord_loadings(mcia_out = mcia_results, omic = "all",
                          absolute = FALSE, descending = TRUE, factor = 1)

test_that("factor filtration", {
    expect_equal(unique(all_pos_1$factor), 1)
})

test_that("all omics being used", {
    expect_equal(sort(as.character(unique(all_pos_1$omic))), # this is a factor
                 c("miRNA", "mrna", "prot"))
})
