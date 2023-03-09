data(NCI60)

# Calculate TV using long method
#mrna_c<-sweep(data_blocks$mrna, 2,colMeans(data_blocks$mrna), "-")
#mirna_c<-sweep(data_blocks$miRNA, 2,colMeans(data_blocks$miRNA), "-")
#prot_c<-sweep(data_blocks$prot, 2,colMeans(data_blocks$prot), "-")

#p1<-prcomp(mrna_c,scale=FALSE)
#p1_s<-p1$sdev
#tv_p1<-sum(p1_s^2)

#p2<-prcomp(mirna_c,scale=FALSE)
#p2_s<-p2$sdev
#tv_p2<-sum(p2_s^2)

#p3<-prcomp(prot_c,scale=FALSE)
#p3_s<-p3$sdev
#tv_p3<-sum(p3_s^2)

#tv_final = list(mrna = tv_p1, mirna= tv_p2, prot= tv_p3)
tv_final<-list(mrna=10152.24,miRNA=1071.079,prot=28843.05)

test_that("get_TV equals standard method", {
  expect_equal(get_TV(data_blocks), tv_final, tolerance=1e-4)
})

