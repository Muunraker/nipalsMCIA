# projection_plot uses base R plots

test_that("projection_plot: checking global figure", {
    expect_no_error(projection_plot(mcia_results, projection = "global",
                                    orders = c(1, 2)))
})

test_that("projection_plot: checking global figure with colors", {
    expect_no_error(projection_plot(mcia_results, projection = "global",
                                    orders = c(1, 2), color_col = "cancerType",
                                    color_pal = meta_colors,
                                    legend_loc = "bottomleft"))
})

test_that("projection_plot: checking all figures", {
    expect_no_error(projection_plot(mcia_results, projection = "all",
                                    orders = c(1, 2)))
})

test_that("projection_plot: checking block figure", {
    expect_no_error(projection_plot(mcia_results, projection = "block",
                                    block_name = "prot",
                                    orders = c(1, 2)))
})

test_that("projection_plot: checking block figure - fail", {
    expect_error(projection_plot(mcia_results, projection = "block",
                                 block_name = "faux-block",
                                 orders = c(1, 2)))
    dev.off()
})
