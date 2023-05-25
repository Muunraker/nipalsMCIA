simple_input<-matrix(c(1,0,0,0,1,0,0,1,0),nrow=3,ncol=3)
temp_df<-simple_input

totsum<-3
colsums<-c(1,1,1)
row_contribs <- c(1,2,0)/totsum

temp_df<-temp_df - row_contribs

temp_df<-t(t(temp_df)*sqrt(colsums/totsum))
block_var<-norm(temp_df,type="F")^2/2
simple_input_cc_preproc<-temp_df/sqrt(block_var)


test_that("cc_preproc works on simple input", {
  expect_equal(simple_input_cc_preproc, cc_preproc(simple_input))
})
