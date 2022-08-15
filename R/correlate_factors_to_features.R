library("ggplot2")
library("reshape2")
library("stringr")
library("ggpubr")
library("gridExtra")
library("dplyr")
library("ComplexHeatmap")
library("circlize")
library('biomaRt')

# Annotation databases used for biological enrichment
path.database <- "data/bio_annotations/c2.cp.reactome.v6.2.symbols.gmt" #REACTOME

# load the factors 
load("results/factors/factorizations.RData")

cytof_fn = 'data/cmi_pb/cytof.pma.day0.proc.tsv'
cytof = read.table(cytof_fn, header=T, row.names = 1)

olink_fn = 'data/cmi_pb/olink.pma.day0.proc.tsv'
olink = read.table(olink_fn, header=T, row.names = 1)

rnaseq_fn = 'data/cmi_pb/rnaseq.pma.day0.proc.tsv'
rnaseq = read.table(rnaseq_fn, header=T, row.names = 1)

features = list(cytof=cytof, olink=olink, rnaseq=rnaseq)

results_folder = "results/factors_v_features/"
dir.create(results_folder, showWarnings = F)

# loading the biomart
ensembl = useMart("ensembl", host="useast.ensembl.org", version =  )
mart <- useDataset("hsapiens_gene_ensembl", ensembl)

# change the metagene ensembl id to hgnc name
for (i in seq(1, length(out$factorizations))){
  genes = rownames(out$factorizations[[i]][[2]][[1]])
  
  print("obtaining hgnc name")
  hgnc_list <- getBM(filters="ensembl_gene_id",
                     attributes=c("ensembl_gene_id", "hgnc_symbol"),
                     values=genes,
                     uniqueRows = F,
                     verbose = 104,
                     mart=mart)
  print("renaming now")
  new_factored_data =  out$factorizations[[i]][[2]][[1]]
  new_factored_data = merge(new_factored_data, hgnc_list, by.x=0, by.y='ensembl_gene_id')
  
  print("running distinct")
  new_factored_data = distinct(new_factored_data, hgnc_symbol, .keep_all = T)
  rownames(new_factored_data) = new_factored_data$hgnc_symbol
  
  print("removing columns")
  new_factored_data = subset(new_factored_data, select = -c(Row.names, hgnc_symbol))
  
  print("saving to a new list in the out variable")
  out$factorizations[[i]][[3]] = list()
  out$factorizations[[i]][[3]][[1]] = new_factored_data
  print("renaming complete")
}
  
############################################
########## Using ggplot2 heatmaps ##########
############################################

# Function: correlates components versus feature data
make_basic_heatmap <- function(feat, feat_data, factor_out){
  
  plot_list = list()
  
  # correlate the tasks and factors
  for (i in seq(1,5)){
    
    # harmonize the data
    method = factor_out$method[[i]]
    factors = factor_out$factorizations[[i]][[1]]
    row.names(factors) = gsub('X', '', row.names(factors))
    shared_samples = intersect(row.names(factors), row.names(feat_data))
    final_factors = factors[shared_samples, ]
    final_features = feat_data[shared_samples, ]
    
    # correlate the factors and task
    factor_corrs = cor(final_factors, final_features)
    factor_corrs = factor_corrs[ , colSums(is.na(factor_corrs)) == 0]
    factor_corrs = melt(factor_corrs)
    colnames(factor_corrs) = c('comp', 'feature', 'comp_val')
    factor_corrs[, 'comp'] = gsub("^.*([0-9]+)$", '\\1', factor_corrs[, 'comp'])
    
    # make a heatmap of the correlations
    p = ggplot(factor_corrs, aes(comp, feature, fill= comp_val)) + 
      geom_tile()
    title = sprintf('%s components versus %s', method, toupper(feat_name))
    p = p + ggtitle(title)
    p = p + theme(plot.title = element_text(size=12)) 
    p = p + scale_fill_distiller(palette = 'RdYlBu', limits=c(-1, 1))
    plot_list[[i]] = p
  }
  return(plot_list)
}

for (feat_name in names(features)){
  
  # generate heatmaps in ggplot
  plot_list = make_basic_heatmap(feat_name, features[[feat_name]], out)

  # save the heatmaps to separate pdf pages
  fn = paste0(results_folder, feat_name, '.pdf')
  pdf(fn, onefile = TRUE)
  for (i in seq(1, length(plot_list))){
    print(plot_list[[i]])
  }
  dev.off()
}

############################################
########## Using ComplexHeatmaps ###########
############################################

# Function: correlates components versus feature data
make_basic_heatmap_v2 <- function(feat, feat_data, factor_out){
  
  plot_list = list()
  
  # correlate the tasks and factors
  for (i in seq(1,5)){
    
    # harmonize the data
    method = factor_out$method[[i]]
    factors = factor_out$factorizations[[i]][[1]]
    row.names(factors) = gsub('X', '', row.names(factors))
    shared_samples = intersect(row.names(factors), row.names(feat_data))
    final_factors = factors[shared_samples, ]
    final_features = feat_data[shared_samples, ]
    
    # correlate the factors and task
    factor_corrs = cor(final_factors, final_features)
    factor_corrs = factor_corrs[ , colSums(is.na(factor_corrs)) == 0]
    rownames(factor_corrs) = gsub("^.*([0-9]+)$", '\\1', rownames(factor_corrs))
    color_func = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    
    # make a heatmap of the correlations
    title = sprintf('%s components versus %s', method, toupper(feat_name))
    show_cols = T
    if (feat_name %in% c('rnaseq', 'olink')){
      show_cols = F
    }
    
    p = Heatmap(factor_corrs, 
            name = "Pearson's R", 
            column_title = title,
            row_title = "Components",
            row_names_gp = gpar(fontsize = 7),
            col = color_func,
            show_column_names = show_cols
    )    
    plot_list[[i]] = p
  }
  return(plot_list)
}

for (feat_name in names(features)){
  
  print(feat_name)
  
  # generate heatmaps in ggplot
  plot_list = make_basic_heatmap_v2(feat_name, features[[feat_name]], out)

  # save the heatmaps to separate pdf pages
  fn = paste0(results_folder, feat_name, '_v3.pdf')
  pdf(fn, onefile = TRUE)
  for (i in seq(1, length(plot_list))){
    print(plot_list[[i]])
  }
  dev.off()
}

############################################
########## Filtering #######################
############################################

# filtering based on the variances of gene expression
# filtering is done using a percentile min and max
rnaseq_variance = rnaseq %>% summarise_all(var)
rnaseq_variance =  t(rnaseq_variance)

# Returns a vector for easy filtering
# filter based on the percentile using pmin and pmax
filter_by_percentile <- function(vec, pmin, pmax){
  var_func = ecdf(rnaseq_variance[, 1])
  bools = list()
  for (i in seq(1, length(vec))){
    p = var_func(vec[i])
    if (p >= pmin & p <= pmax){
      bools[[i]] = T
    }
    else{
      bools[[i]] = F
    }
  }
  bools = unlist(bools)
  return(bools)
}
flt_bools = filter_by_percentile(rnaseq_variance[, 1], 0.995, 1)
flt_genes = rownames(rnaseq_variance)[flt_bools]
flt_rnaseq = rnaseq[,flt_genes]

# generate heatmaps in ggplot
plot_list = make_basic_heatmap_v2('rnaseq', flt_rnaseq, out)

# save the heatmaps to separate pdf pages
fn = paste0(results_folder, 'rnaseq_flt', '_v3.pdf')
pdf(fn, onefile = TRUE)
for (i in seq(1, length(plot_list))){
  print(plot_list[[i]])
}
dev.off()


########################################
######## biological comparison #########
########################################

## Perform biological annotation-based comparison 

library("fgsea", quietly = TRUE)

## Perform biological annotation-based comparison 
## INPUTS:
# factorizations = already computed factorizations
# path.database = path to a GMT annotation file
# pval.thr = p-value threshold (default to 0.05)
## OUPUTS: a list containing output values
# selectivity = Selectivity (fraction of significant annotations per all significant factors)
# nonZeroFacs = Number of unique factors that were significant at least once
# total_pathways = Number of clinical annotations with at least one significant factor component
biological_comparison <- function(factorizations, path.database, pval.thr=0.05){
    
    # Load annotation database
    pathways <- gmtPathways(path.database)
    
    # Containers to report results
    report_number <- numeric(0)
    report_nnzero <- numeric(0)
    report_select <- numeric(0)
    
    # For each factorization method
    for(i in 1:length(factorizations)){
        
        cat(paste0("Studying factor: ", i, "\n"))
        
        # Extract metagenes found by factorization method
        cat("# Extract metagenes found by factorization method\n")
        metagenes <- factorizations[[i]][[3]][[1]]
        
        #cat("metagenes 1:\n")
        #print(metagenes[1:10, 1:10])
        
        # Number of factors
        num.factors <- ncol(metagenes)
        cat(paste0("# Number of factors: ", num.factors, "\n"))
        
        #cat("metagenes 2:\n")
        #print(metagenes[1:10, 1:10])
        
        # Rename columns
        cat("# Rename columns\n")
        colnames(metagenes) <- 1:num.factors
        
        # Rename rows to remove "|" characters and keep only the gene name before
        cat("# Rename rows to remove '|' characters and keep only the gene name before\n")
        rownames(metagenes) <- gsub("\\|",".",rownames(metagenes))
        rownames(metagenes) <- gsub("\\..*","",rownames(metagenes))
        
        # Remove duplicated gene names that could confuse fgsea
        cat("# Remove duplicated gene names that could confuse fgsea\n")
        duplicated_names <- unique(rownames(metagenes)[duplicated(rownames(metagenes))])
        metagenes <- metagenes[!(rownames(metagenes) %in% duplicated_names), ]
        
        #cat("metagenes 3:\n")
        #print(metagenes[1:10, 1:10])
        
        # Variables
        min_pval <- numeric(0)
        path <- numeric(0)
        n <- 0
        
        # Calculate biological annotation enrichment.
        cat("# Calculate biological annotation enrichment.\n")
        # For each factor,
        for(j in 1:num.factors){
          
            # Assign gene names
            cat("\t# Assign gene names\n")
            rnk <- setNames(as.matrix(metagenes[,j]), rownames(metagenes))
            
            # Compute fgsea
            cat("\t# Compute fgsea\n")
            #fgseaRes <- fgsea(pathways, rnk, minSize=15, maxSize=500, nperm=1000)
            fgseaRes <- fgseaMultilevel(pathways, rnk, nPermSimple = 10000,minSize=15, maxSize=500)
            
            # If at least one pathway is significant
            cat("\t# If at least one pathway is significant\n")
            
            print(paste0("fgseaRes: ", fgseaRes))
            
            
            if(sum(fgseaRes$padj < pval.thr)!=0){
                
                # Count this factor
                n <- n+1
                
                # Keep min adjusted p-value
                min_pval <- rbind(min_pval, min(fgseaRes$padj))
                
                # Keep names of significant pathways
                path <- c(path, fgseaRes[fgseaRes$padj<pval.thr, "pathway"])
                
            } else {
                min_pval <- rbind(min_pval, NA)
            }
        }

        # Report number of unique significant pathways  
        cat("# Report number of unique significant pathways\n")
        if(length(path)==0){
            report_number <- rbind(report_number, 0)
        }else{
            report_number <- rbind(report_number, length(unique(path)))
        }
        
        # Report selectivity 
        cat("# Report selectivity\n")
        if(length(unique(path))==0){
            report_select <- rbind(report_select, NA)
        }else{
            al<-length(unique(path))/length(path)
            fl<-length(which(!is.na(min_pval)))/length(path)
            report_select <- rbind(report_select, (al+fl)/2)
        }
        
        # Report number of factors associated with at least one significant pathway
        cat("# Report number of factors associated with at least one significant pathway\n")
        report_nnzero <- rbind(report_nnzero, n)    
    
        cat("\n")
    }
    
    # reporting the selectivity, number of factors associated with at 
    # least one significant pathway, and the total number of significant
    # pathways
    out <- data.frame(selectivity=report_select,
                      nonZeroFacs=report_nnzero,
                      total_pathways=report_number)
    return(out)
}

bio_comp = biological_comparison(out$factorizations, path.database, pval.thr=0.05)
fn = paste0(results_folder, 'report.tsv')
write.table(bio_comp, file=fn, quote = F, sep="\t")










