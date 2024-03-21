test_that("figure generation", {
    #dev.new()
    #par(mar = c(2, 2, 2, 2))
    expect_no_error(vis_load_plot(mcia_results, axes = c(1, 2)))
    #dev.off()
})

test_that("incorrect axes", {
    #dev.new()
    #par(mar = c(2, 2, 2, 2))
    expect_error(vis_load_plot(mcia_results, axes = c(1, 10)))
    #dev.off()
})

test_that("default axes", {
    expect_no_error(vis_load_plot(mcia_results))
})
