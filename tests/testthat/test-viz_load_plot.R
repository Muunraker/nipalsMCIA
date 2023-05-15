test_that("viz_load_plot: checking figure", {
    expect_no_error(vis_load_plot(mcia_results,
                                  axes = c(1, 2),
                                  colors_omics = colors_omics)
    )
})

test_that("viz_load_plot: checking incorrect axes", {
    expect_error(vis_load_plot(mcia_results,
                                  axes = c(1, 10),
                                  colors_omics = colors_omics)
    )
})