data("NCI60") # loads in "data_blocks" which has 3 types of omics and metadata
mcia_results <- nipals_multiblock(data_blocks, metadata = metadata_NCI60,
                                  preproc_method = "colprofile",
                                  num_PCs = 2, tol = 1e-12, plots = "none")

# TO DO: test that the colors are actually being set properly (not just length)

test_that("generate omics color palette", {
  expect_length(get_colors(mcia_results), 3)
})

test_that("use given omics color palette", {
    expect_length(get_colors(mcia_results,
                             color_pal = c("red", "green", "blue")), 3)
})

test_that("generate metadata color palette", {
    expect_length(get_metadata_colors(mcia_results), 3)
})
