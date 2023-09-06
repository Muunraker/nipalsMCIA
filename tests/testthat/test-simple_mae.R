sce1=matrix(rpois(100, lambda = 10), ncol = 10)
sce2=matrix(rpois(100, lambda = 10), ncol = 10)
sce_list<-list(sce1,sce2)

test_that("All omics must have sample names", 
          {expect_error(simple_mae(sce_list))
})