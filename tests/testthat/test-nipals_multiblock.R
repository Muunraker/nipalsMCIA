# test dimension of global loadings (two omics with equal features)
test_that("Number of block score weights correct (random)", {
  # test data with 3 samples, 2 omics, 5 features per omic
  test_data <- list(omic1 = data.frame((rand(3, 5))),
                    omic2 = data.frame((rand(3, 5))))
  test <- nipals_multiblock(test_data, plots = "none",
                            num_PCs = 4, tol = 1e-9)

  expect_equal(c(10, 4), dim(test$global_loadings))
})

test_that("Number of block score weights correct (NCI60)", {
  expect_equal(c(20448, 2), dim(mcia_results$global_loadings))
})

# testing magnitude of top three eigenvalues of NCI60 dataset
test_that("Top eigenvalues correct", {
  expect_equal(all(mcia_results3$eigvals > 0.21), TRUE)
})

# expect error if number of samples is wrong
test_that("Mismatched samples should throw error", {
  test_data <- list(omic1 = data.frame((rand(2, 5))),
                    omic2 = data.frame((rand(3, 5))))

  expect_error(nipals_multiblock(test_data, plots = "none",
                                 num_PCs = 4, tol = 1e-9))
})
