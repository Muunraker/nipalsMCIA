# extract mRNA global loadings
mrna_gfscores <- mcia_results$global_loadings
mrna_rows <- str_detect(row.names(mrna_gfscores), "_mrna")
mrna_gfscores <- mrna_gfscores[mrna_rows, ]

# rename rows to contain HUGO based gene symbols
row.names(mrna_gfscores) <- str_remove(rownames(mrna_gfscores), "_[0-9]*_.*")

# load pathway data
path.database <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/6.2/c2.cp.reactome.v6.2.symbols.gmt"
pathways <- fgsea::gmtPathways(gmt.file = path.database)

test_that("gsea_report", {
    # generate the GSEA report
    expect_no_error(geneset_report <- gsea_report(metagenes = mrna_gfscores,
                                                  path.database,
                                                  factors = seq(1),
                                                  pval.thr = 0.05,
                                                  nproc = 1))
})
