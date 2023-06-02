bw_heatmap <- block_weights_heatmap(mcia_results)
bs_weights <- as.matrix(data.frame(mcia_results$block_score_weights))
colnames(bs_weights) <- seq_len(ncol(bs_weights))

test_that("block_weights_heatmap: checking figure", {
    expect_no_error(block_weights_heatmap(mcia_results))
})

test_that("input is correct", {
  expect_equal(bw_heatmap@matrix, bs_weights)
})
