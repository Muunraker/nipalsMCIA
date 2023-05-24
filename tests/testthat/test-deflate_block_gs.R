test_df <- cbind(c(0, 0, 1, 1, 0), c(1, 0, 0, 0, 0), c(0, 0, 0, 1, 1))
global_score <- c(0, 0, 1, 0, 0)

normed_gs <-  global_score/norm(global_score, type = "2")
deflate_gs <- test_df - tcrossprod(normed_gs) %*% as.matrix(test_df)

test_that("test output",  {
  expect_equal(deflate_gs, deflate_block_gs(test_df, global_score))
})

