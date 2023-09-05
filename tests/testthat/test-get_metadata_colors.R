# uses the default viridis palette
test_that("generate metadata color palette", {
  expect_length(get_metadata_colors(mcia_results), 3)
})

# tries a different type of palette
test_that("check metadata color palette", {
  expect_equal(
    get_metadata_colors(mcia_results, color_pal = scales::brewer_pal,
                        color_pal_params = list()),
    c("CNS" = "#DEEBF7", "Leukemia" = "#9ECAE1", "Melanoma" = "#3182BD")
  )
})

# tries a different type of palette with a different option
test_that("check metadata color palette", {
  expect_equal(
    get_metadata_colors(mcia_results, color_pal = scales::brewer_pal,
                        color_pal_params = list(palette = "Greens",
                                                direction = -1)),
    c("CNS" = "#31A354", "Leukemia" = "#A1D99B", "Melanoma" = "#E5F5E0")
  )
})
