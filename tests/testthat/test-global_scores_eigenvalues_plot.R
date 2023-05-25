# global_scores_eigenvalues_plot uses base R plots
test_that("checking figure", {
    dev.new()
    par(mar = c(2, 2, 2, 2))
    expect_no_error(global_scores_eigenvalues_plot(mcia_results))
    dev.off()
})
