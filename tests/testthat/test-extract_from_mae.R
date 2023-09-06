test_that("The nipalsMCIA algorithm is designed for analysis of > 1 omics experiments", 
          {expect_error(extract_from_mae(data_blocks_mae, subset="prot"))
})

test_that("Please provide appropriate subset_data list", 
          {expect_error(extract_from_mae(data_blocks_mae, subset="chrom"))
          })