test_that("block preprocessing by number of columns", {
    prot_num_cols <- prot / ncol(prot)

    expect_equal(block_preproc(prot, "num_cols"), prot_num_cols)
})


test_that("block preprocessing by unit variance", {
    block_var <- norm(as.matrix(prot), type = "F")^2 / (max(1, nrow(prot) - 1))
    prot_unit_var <- as.matrix(prot) * (1 / sqrt(block_var))
    
    expect_equal(block_preproc(as.matrix(prot),"unit_var"),prot_unit_var)
})

test_that("block preprocessing by largest_sv", {
  svdres <- RSpectra::svds(as.matrix(prot),1)
  prot_svd <- as.matrix(prot)*(1/svdres$d)
  
  expect_equal(block_preproc(as.matrix(prot),"largest_sv"),prot_svd)
})

test_that("block preprocessing by none", {
  expect_equal(block_preproc(prot,"none"),prot)
})