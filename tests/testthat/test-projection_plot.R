test_that("projection_plot: checking global figure", {
    dev.new()
    expect_no_error(projection_plot(mcia_results,
                    projection = "global",
                    orders = c(1, 2)))
    dev.off()

})