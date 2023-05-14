gs_heatmap<-global_scores_heatmap(mcia_results)

test_that("Name is correct", {
  expect_equal(gs_heatmap@name,"GS Score")
})
