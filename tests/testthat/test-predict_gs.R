preds <- predict_gs(mcia_results, data_blocks_mae)

test_that("correct dimensions", {
  expect_equal(nrow(preds), nrow(mcia_results@metadata))
})
