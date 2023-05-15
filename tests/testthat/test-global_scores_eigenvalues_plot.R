test_that("global_scores_heatmap: checking figure", {
    dev.new()
    expect_no_error(global_scores_eigenvalues_plot(mcia_results))
    dev.off()
})
