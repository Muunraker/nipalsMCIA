test_that("checking figure", {
    expect_no_error(global_scores_heatmap(mcia_results,
                                         color_col="cancerType")
    )
})

test_that("checking figure with custom color palette", {
    expect_no_error(global_scores_heatmap(mcia_results,
                                          color_col="cancerType",
                                          color_pal = meta_colors)
    )
})


test_that("global_scores_heatmap: color column not in metadata", {
    expect_error(global_scores_heatmap(mcia_results_no_meta,
                                         color_col="faux-column")
    )
})

gs_heatmap<-global_scores_heatmap(mcia_results)

test_that("Name is correct", {
  expect_equal(gs_heatmap@name,"GS Score")
})
