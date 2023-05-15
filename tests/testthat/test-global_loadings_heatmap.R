# test_that("metadata not present", {
#     expect_error(global_loadings_heatmap(mcia_results_no_meta,
#                                          color_col="faux-column")
#     )
# })
# 
# test_that("omic not present in blocks", {
#     expect_error(global_loadings_heatmap(mcia_results,
#                                          data_blocks,
#                                          "faux-omic")
#     )
# })
# 
# test_that("checking figure creation", {
#     expect_no_error(global_loadings_heatmap(mcia_results,
#                                            color_col="cancerType")
#     )
# })
