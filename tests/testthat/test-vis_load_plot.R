test_that("figure generation", {
    #dev.new()
    #par(mar = c(2, 2, 2, 2))
    expect_no_error(vis_load_plot(mcia_results, axes = c(1, 2)))
    #dev.off()
})

test_that("incorrect axes", {
    #dev.new()
    #par(mar = c(2, 2, 2, 2))
    expect_error(vis_load_plot(mcia_results, axes = c(1, 10)))
    #dev.off()
})

# when we add in a default color palette to this function, we should switch this to expect_no_error
test_that("color scheme", {
    expect_error(vis_load_plot(mcia_out = mcia_results, axes = c(1, 2)))
})

test_that("default axes", {
    expect_no_error(vis_load_plot(mcia_results))
})
