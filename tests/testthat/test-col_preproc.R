data(NCI60)

# Extract just proteomics dataset
prot = data_blocks$prot

# Standardized proteomics dataset
prot = as.matrix(prot)
prot_stand = scale(prot)

test_that("standardized column preprocessing works", {
  expect_equal(col_preproc(prot,"standardized"), prot_stand)
})
