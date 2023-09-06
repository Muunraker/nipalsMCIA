# "data_blocks" (loaded with `data(NCI60)` during setup)
# has 3 types of omics and metadata

test_that("generate omics color palette", {
    expect_length(get_colors(mcia_results), 3)
})

# uses the default viridis palette
test_that("check omics color palette 1", {
    expect_equal(
        get_colors(mcia_results), # option = "D"
        c("mrna" = "#440154FF", "miRNA" = "#21908CFF", "prot" = "#FDE725FF")
    )
})

# tries a different type of viridis palette
test_that("check omics color palette 2", {
    expect_equal(
        get_colors(mcia_results, color_pal_params = list(option = "A")),
        c("mrna" = "#000004FF", "miRNA" = "#B63679FF", "prot" = "#FCFDBFFF")
    )
})

# tries a different type of palette
test_that("check omics color palette 3", {
    expect_equal(
        get_colors(mcia_results, color_pal = scales::hue_pal),
        c("mrna" = "#F8766D", "miRNA" = "#00BA38", "prot" = "#619CFF")
    )
})

# tries a different type of palette with a different option
test_that("check omics color palette 4", {
    expect_equal(
        get_colors(mcia_results, color_pal = scales::hue_pal,
                   color_pal_params = list(h = c(0, 90))),
        c("mrna" = "#FF6C91", "miRNA" = "#DE8C00", "prot" = "#9DA700")
    )
})

test_that("use given omics color palette", {
    expect_length(get_colors(mcia_results,
                             color_pal = c("red", "green", "blue")), 3)
})


