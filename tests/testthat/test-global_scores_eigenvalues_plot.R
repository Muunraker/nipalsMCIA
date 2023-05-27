# global_scores_eigenvalues_plot uses base R plots
test_that("checking figure", {
    expect_no_error(global_scores_eigenvalues_plot(mcia_results))
})
