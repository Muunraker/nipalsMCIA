# projection_plot uses base R plots

test_that("projection_plot: checking global figure", {
    dev.new()
    par(mar = c(2, 2, 2, 2)) 
    expect_no_error(projection_plot(mcia_results,
                    projection = "global",
                    orders = c(1, 2)))
    dev.off()
})

test_that("projection_plot: checking global figure with colors", {
    dev.new()
    par(mar = c(2, 2, 2, 2)) 
    expect_no_error(projection_plot(mcia_results, 
                                    projection = "global", 
                                    orders = c(1, 2),
                                    color_col = "cancerType",
                                    color_pal = meta_colors,
                                    legend_loc = "bottomleft")
    )
    dev.off()
})

test_that("projection_plot: checking all figure", {
    dev.new()
    par(mar = c(2, 2, 2, 2)) 
    expect_no_error(projection_plot(mcia_results,
                                    projection = "all",
                                    orders = c(1, 2)))
    dev.off()
})

test_that("projection_plot: checking block figure", {
    dev.new()
    par(mar = c(2, 2, 2, 2)) 
    expect_no_error(projection_plot(mcia_results,
                                    projection = "block",
                                    block_name = "prot",
                                    orders = c(1, 2)))
    dev.off()
})

test_that("projection_plot: checking block figure - fail", {
    dev.new()
    par(mar = c(2, 2, 2, 2)) 
    expect_error(projection_plot(mcia_results,
                                    projection = "block",
                                    block_name = "faux-block",
                                    orders = c(1, 2)))
    dev.off()
})