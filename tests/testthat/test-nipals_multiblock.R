# test if # PCs generated matches # PCs specified
test_that("Number of block score weights correct", {
  test_data<-list(omic1=data.frame((rand(3,5))),omic2=data.frame((rand(3,5))))
  test <- nipals_multiblock(test_data,plots = 'none',num_PCs = 10, tol=1e-9)
  expect_equal(2,dim(test$block_score_weights)[1])
})

# test if dimensions of block/global loadings are correct

# test if works in edge case: # samples = # variables

# test if output works for # samples > # variables

# test if output expected for one block(?) 


