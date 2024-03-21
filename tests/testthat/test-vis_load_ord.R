mrna_pos_1 <- ord_loadings(mcia_results, omic = "mrna",
                           absolute = FALSE, descending = TRUE, factor = 1)
mrna_pos_1_vis <- vis_load_ord(mcia_results, omic = "mrna")

test_that("omic is correct", {
  expect_equal(as.character(unique(mrna_pos_1_vis$data$omic)[1]), "mrna")
})

test_that("factor is correct", {
  expect_equal(as.character(unique(mrna_pos_1_vis$data$factor)[1]), "1")
})

test_that("no errors produced", {
  expect_no_error(vis_load_ord(mcia_results, omic = "mrna"))
})
