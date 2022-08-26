#library("data.table", quietly = T)
#library("ggplot2", quietly=T)
#library("reshape2", quietly=T)
#library("stringr", quietly=T)
#library("ggpubr", quietly=T)
#library("gridExtra", quietly=T)
#library("dplyr", quietly=T)
#library("ComplexHeatmap", quietly=T)
#library("circlize", quietly=T)
#library('biomaRt', quietly=T)
#library("fgsea", quietly = T)
#
## load the factors 
#load("results/factors/factorizations.RData")
#
## set the results folder; needs to come after the previous load
## or else it is overridden
#results_folder = "results/pathway_analysis/"
#dir.create(results_folder, showWarnings = F)
#
## Annotation databases used for biological enrichment
#path.database <- "data/bio_annotations/c2.cp.reactome.v6.2.symbols.gmt" #REACTOME
#
## loading the biomart
#ensembl = useMart("ensembl", host="useast.ensembl.org", version =  )
#mart <- useDataset("hsapiens_gene_ensembl", ensembl)
#
## change the metagene ensembl id to hgnc name
#for (i in seq(1, length(out$factorizations))){
#  genes = rownames(out$factorizations[[i]][[2]][[1]])
#  
#  print("obtaining hgnc name")
#  hgnc_list <- getBM(filters="ensembl_gene_id",
#                     attributes=c("ensembl_gene_id", "hgnc_symbol"),
#                     values=genes,
#                     uniqueRows = F,
#                     verbose = 104,
#                     mart=mart)
#  print("renaming now")
#  new_factored_data =  out$factorizations[[i]][[2]][[1]]
#  new_factored_data = merge(new_factored_data, hgnc_list, by.x=0, by.y='ensembl_gene_id')
#  
#  print("running distinct")
#  new_factored_data = distinct(new_factored_data, hgnc_symbol, .keep_all = T)
#  rownames(new_factored_data) = new_factored_data$hgnc_symbol
#  
#  print("removing columns")
#  new_factored_data = subset(new_factored_data, select = -c(Row.names, hgnc_symbol))
#  
#  print("saving to a new list in the out variable")
#  out$factorizations[[i]][[3]] = list()
#  out$factorizations[[i]][[3]][[1]] = new_factored_data
#  print("renaming complete")
#}
#
### Perform biological annotation-based comparison 
### INPUTS:
## factorizations = already computed factorizations
## path.database = path to a GMT annotation file
## pval.thr = p-value threshold (default to 0.05)
### OUPUTS: a list containing output values
## selectivity = Selectivity (fraction of significant annotations per all significant factors)
## nonZeroFacs = Number of unique factors that were significant at least once
## total_pathways = Number of clinical annotations with at least one significant factor component
#biological_comparison <- function(factorizations, path.database, pval.thr=0.05){
#    
#    # Load annotation database
#    pathways <- gmtPathways(path.database)
#    
#    # Containers to report results
#    report_number <- numeric(0)
#    report_nnzero <- numeric(0)
#    report_select <- numeric(0)
#    
#    # For each factorization method
#    for(i in 1:length(factorizations)){
#        
#        cat(paste0("Studying factor: ", i, "\n"))
#        
#        # Extract metagenes found by factorization method
#        cat("# Extract metagenes found by factorization method\n")
#        metagenes <- factorizations[[i]][[3]][[1]]
#        
#        # Number of factors
#        num.factors <- ncol(metagenes)
#        cat(paste0("# Number of factors: ", num.factors, "\n"))
#        
#        # Rename columns
#        cat("# Rename columns\n")
#        colnames(metagenes) <- 1:num.factors
#        
#        # Rename rows to remove "|" characters and keep only the gene name before
#        cat("# Rename rows to remove '|' characters and keep only the gene name before\n")
#        rownames(metagenes) <- gsub("\\|",".",rownames(metagenes))
#        rownames(metagenes) <- gsub("\\..*","",rownames(metagenes))
#        
#        # Remove duplicated gene names that could confuse fgsea
#        cat("# Remove duplicated gene names that could confuse fgsea\n")
#        duplicated_names <- unique(rownames(metagenes)[duplicated(rownames(metagenes))])
#        metagenes <- metagenes[!(rownames(metagenes) %in% duplicated_names), ]
#        
#        # Variables
#        min_pval <- numeric(0)
#        path <- numeric(0)
#        n <- 0
#        
#        # Calculate biological annotation enrichment.
#        cat("# Calculate biological annotation enrichment.\n")
#        # For each factor,
#        for(j in 1:num.factors){
#          
#            # Assign gene names
#            cat("\t# Assign gene names\n")
#            rnk <- setNames(as.matrix(metagenes[,j]), rownames(metagenes))
#            
#            # Compute fgsea
#            cat("\t# Compute fgsea\n")
#            fgseaRes <- fgsea(pathways, rnk, minSize=15, maxSize=500)
#            
#            # If at least one pathway is significant
#            cat("\t# If at least one pathway is significant\n")
#            
#            print(paste("colnames(fgseaRes):", colnames(fgseaRes)))
#            print(paste0("fgseaRes: ", fgseaRes[1:10, 1:7]))
#            print(paste0("fgseaRes$padj: ", fgseaRes$padj))
#            
#            # count this factor if in some significant pathway
#            pval_check = sum(fgseaRes$padj < pval.thr, na.rm=T)
#            if(pval_check != 0){
#                
#                # Count this factor
#                n <- n + 1
#                
#                # Keep min adjusted p-value
#                min_pval <- rbind(min_pval, min(fgseaRes$padj))
#                
#                # Keep names of significant pathways
#                path <- c(path, fgseaRes[fgseaRes$padj<pval.thr, "pathway"])
#                
#            } else {
#                min_pval <- rbind(min_pval, NA)
#            }
#        }
#
#        # Report number of unique significant pathways  
#        cat("# Report number of unique significant pathways\n")
#        if(length(path)==0){
#            report_number <- rbind(report_number, 0)
#        }
#        else{
#            report_number <- rbind(report_number, length(unique(path)))
#            fwrite(list(path), file = paste0(results_folder, out$method[[i]], "_paths.txt"))
#        }
#        
#        # Report selectivity 
#        cat("# Report selectivity\n")
#        if(length(unique(path))==0){
#            report_select <- rbind(report_select, NA)
#        }
#        else{
#            al<-length(unique(path))/length(path)
#            fl<-length(which(!is.na(min_pval)))/length(path)
#            report_select <- rbind(report_select, (al+fl)/2)
#        }
#        
#        # Report number of factors associated with at least one significant pathway
#        cat("# Report number of factors associated with at least one significant pathway\n")
#        report_nnzero <- rbind(report_nnzero, n)    
#        cat("\n")
#    }
#    
#    # reporting the selectivity, number of factors associated with at 
#    # least one significant pathway, and the total number of significant
#    # pathways
#    report <- data.frame(selectivity=report_select,
#                      nonZeroFacs=report_nnzero,
#                      total_pathways=report_number)
#    return(report)
#}
#
#bio_comp = biological_comparison(out$factorizations, path.database, pval.thr=0.05)
#rownames(bio_comp) = out$method
#fn = paste0(results_folder, 'report.tsv')
#write.table(bio_comp, file=fn, quote = F, sep="\t")
