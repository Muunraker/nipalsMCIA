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

#' General purpose function for plotting a correlation heatmap of factors versus other
#'
#' @description Plots a heatmap of MCIA factos versus other vectors such as
#' the original feature vectors (membership analysis) or versus tasks to 
#' investigate the predictive power of these factors to some other vector of data.
#' @param block_matrix the list of matrices used to calculate the MCIA factors
#' @param mcia_factors the MCIA factor object calculated from block_matrix
#' @return the ggplot2 object
#' @export
plot_factors_v_features <- function(feat, feat_data, factor_out){
  
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

#' General purpose function for plotting a correlation heatmap of factors versus other
#'
#' @description Plots a heatmap of MCIA factos versus other vectors such as
#' the original feature vectors (membership analysis) or versus tasks to 
#' investigate the predictive power of these factors to some other vector of data.
#' @param block_matrix the list of matrices used to calculate the MCIA factors
#' @param mcia_factors the MCIA factor object calculated from block_matrix
#' @return the ggplot2 object
#' @export
plot_factors_v_features_v2 <- function(feat, feat_data, factor_out){

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

#' Perform biological annotation-based comparison 
#'
#' @description Runs fgsea for the input gene list
#' @param genes a vector of gene names according to HUGO nomenclature 
#' @param factorizations = already computed factorizations
#' @param path.database = path to a GMT annotation file
#' @param pval.thr = p-value threshold (default to 0.05)
#' @return selectivity Selectivity (fraction of significant annotations per all significant factors)
#' @return nonZeroFacs Number of unique factors that were significant at least once
#' @return total_pathways = Number of clinical annotations with at least one significant factor component
#' @export
biological_comparison <- function(factorizations, path.database, pval.thr=0.05){

    library("fgsea", quietly = TRUE)
    
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

#bio_comp = biological_comparison(out$factorizations, path.database, pval.thr=0.05)
#fn = paste0(results_folder, 'report.tsv')
#write.table(bio_comp, file=fn, quote = F, sep="\t")










