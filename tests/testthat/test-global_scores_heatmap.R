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


test_that("color column not in metadata", {
    expect_error(global_scores_heatmap(mcia_results_no_meta,
                                         color_col="faux-column")
    )
})


test_that("checking the ggplot object with inputs", {

    # capture the plotting function
    p <- global_scores_heatmap(mcia_results, color_col="cancerType")
    expect_equal(matrix(p@matrix), matrix(mcia_results$global_scores))
    expect_equal("cancerType", p@right_annotation@anno_list$ColorType@label)

})


gs_heatmap <- global_scores_heatmap(mcia_results)
test_that("Name is correct", {
  expect_equal(gs_heatmap@name,"GS Score")
})
