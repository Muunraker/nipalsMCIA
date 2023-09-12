# test dimension of global loadings (two omics with equal features)
test_that("Number of block score weights correct (random)", {
  # test data with 3 samples, 2 omics, 5 features per omic
  omic1 <- data.frame((rand(3, 5)))
  omic2 <- data.frame((rand(3, 5)))
  rownames(omic1) <- rownames(omic2) <- c("s1", "s2", "s3")
  colnames(omic1) <- c("f11", "f21", "f31", "f41", "f51")
  colnames(omic2) <- c("f12", "f22", "f32", "f42", "f52")

  test_data <- list("omic1" = omic1, "omic2" = omic2)
  test_data_mae <- simple_mae(test_data, row_format = "sample")

  test <- nipals_multiblock(test_data_mae, plots = "none",
                            num_PCs = 4, tol = 1e-9)

  expect_equal(c(10, 4), dim(test@global_loadings))
})

test_that("Number of block score weights correct (NCI60)", {
  expect_equal(c(20448, 2), dim(mcia_results@global_loadings))
})

# testing magnitude of top three eigenvalues of NCI60 dataset
test_that("Top eigenvalues correct", {
  expect_equal(all(mcia_results3@eigvals > 0.21), TRUE)
})

# expect error if number of samples is wrong
test_that("Mismatched samples should throw error", {
  test_data <- list(omic1 = data.frame((rand(2, 5))),
                    omic2 = data.frame((rand(3, 5))))

  expect_error(nipals_multiblock(test_data, plots = "none",
                                 num_PCs = 4, tol = 1e-9))
})
