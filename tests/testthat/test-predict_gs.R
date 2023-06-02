preds <- predict_gs(mcia_results, data_blocks)

test_that("correct dimensions", {
  expect_equal(nrow(preds), nrow(mcia_results$metadata))
})
