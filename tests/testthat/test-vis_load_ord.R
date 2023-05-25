mrna_pos_1 <- ord_loadings(mcia_out = mcia_results, omic = "mrna",
                           absolute = FALSE, descending = TRUE, factor = 1)
mrna_pos_1_vis <- vis_load_ord(gl_f_ord = mrna_pos_1, omic_name = "mrna",
                               colors_omics = colors_omics)

test_that("omic is correct", {
  expect_equal(as.character(unique(mrna_pos_1_vis$data$omic)[1]), "mrna")
})

test_that("factor is correct", {
  expect_equal(as.character(unique(mrna_pos_1_vis$data$factor)[1]), "1")
})

test_that("no errors produced", {
  dev.new()
  par(mar = c(2, 2, 2, 2))
  expect_no_error(vis_load_ord(gl_f_ord = mrna_pos_1, omic_name = "mrna",
                               colors_omics = colors_omics))
  dev.off()
})
