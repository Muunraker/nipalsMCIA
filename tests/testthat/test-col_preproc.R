test_that("standardized column preprocessing", {
    prot_stand <- scale(as.matrix(prot), center = TRUE, scale = TRUE)

    expect_equal(col_preproc(prot, "standardized"), prot_stand)
})

test_that("centered only column preprocessing", {
    prot_centered <- scale(as.matrix(prot), center = TRUE, scale = FALSE)

    expect_equal(col_preproc(prot, "centered_only"), prot_centered)
})
