<<<<<<< HEAD
library("ggplot2")
library("reshape2")
library("ComplexHeatmap")
library("circlize")

#library("dplyr")
#library('biomaRt')

#' Plotting a heatmap of global_loadings
#'
#' @description Plots a heatmap of MCIA global_loadings
#' @param global_loadings the global_loadings matrix after running MCIA
#' @param data_blocks the list of matrices used to calculate the MCIA factors
#' @param omic_name the name of the omic that should be plot, should be in data_blocks
#' @param select_features a vector of numbers to filter features 
#' @return plot a plot object from ComplexHeatmap
#' @export
global_loadings_heatmap_ComplexHeatmap <- function(global_loadings,
                                                   data_blocks,
                                                   omic_name, 
                                                   select_factors=NULL,
                                                   select_features=NULL){
    
        # get the index of the omics we need
        omics_index = which(names(data_blocks) == omic_name)
        
        # find at what index we should start in the global loadings
        if (omics_index == 1){
            global_start = 1 
        } else {
            i = 1
            global_start = 0
            while (i != omics_index){
                global_start = global_start + ncol(data_blocks[[i]])
                i = i + 1
            }
        }
        
        # global_end is global_start plus the number or rows in the 
        # data_blocks of omic_name
        global_end = global_start + ncol(data_blocks[[omics_index]]) - 1
        
        # extract data for current omic
        # transpose so rows are the factors and columns are the features
        omic_data = t(global_loadings[seq(global_start, global_end),])
        
        # rename rows and columns for plotting
        rownames(omic_data) = paste0("Factor", seq(1, nrow(omic_data)))
        colnames(omic_data) = colnames(data_blocks[[omic_name]])
        
        # make a heatmap of the correlations
        #color_func = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
        coltitle = sprintf(sprintf("%s Features", omic_name))
        
        # filter factors based on select_factors 
        if (!is.null(select_factors)){
            omic_data = omic_data[select_factors,]
        }
        
        # filter features based on select_features 
        if (!is.null(select_features)){
            omic_data = omic_data[, select_features]
        }
        
        # plot the heatmap
        plot = Heatmap(omic_data, 
                name = "Gobal Loading", 
                column_title = coltitle,
                row_title = "Factors",
                row_names_gp = gpar(fontsize = 7),
                show_column_names = T,
                show_row_names = T,
                row_names_side = "right",
                show_column_dend = F,
                cluster_columns = F,
                cluster_rows = F
        )    
      return(plot)
}


#' Plot a correlation heatmap of global scores versus features using ComplexHeatmap
#'
#' @description Plots a heatmap of MCIA factos versus other vectors such as
#' the original feature vectors (membership analysis) or versus tasks to 
#' investigate the predictive power of these factors to some other vector of data.
#' @param global_scores the MCIA factor object calculated from block_matrix
#' @param feature_mat the input matrix used to calculate the MCIA factors
#' @param feature_name the name of the input
#' @param select_features a vector of numbers to filter features 
#' @return plot a plot object from ComplexHeatmap
#' @export
corr_heatmap_fvl_ComplexHeatmap <- function(global_scores,
                                            feature_mat,
                                            feature_name="",
                                            select_features=NULL){
        
        # correlate the factors and task
        lf_corrs = cor(global_scores, feature_mat)
        rownames(lf_corrs) = paste("Factor", seq(1, nrow(lf_corrs)))
        
        # filter features basd on select_features 
        if (!is.null(select_features)){
            lf_corrs = lf_corrs[, select_features]
        }
        
        # make a heatmap of the correlations
        color_func = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
        title = sprintf('%s Features', feature_name)
        plot = Heatmap(lf_corrs, 
                name = "Pearson's R", 
                column_title = title,
                row_title = "Factors",
                row_names_gp = gpar(fontsize = 7),
                col = color_func,
                show_column_names = T,
                show_row_names = T,
                row_names_side = "right",
                cluster_rows = F,
                show_column_dend = F,
                cluster_columns = F
        )    
      return(plot)
}


#' Plot a correlation heatmap of global scores versus features using ggplot2
#'
#' @description Plots a heatmap of MCIA factos versus other vectors such as
#' the original feature vectors (membership analysis) or versus tasks to 
#' investigate the predictive power of these factors to some other vector of data.
#' @param global_scores the MCIA factor object calculated from block_matrix
#' @param feature_mat the input matrix used to calculate the MCIA factors
#' @param feature_name the name of the input
#' @param select_features a vector of numbers to filter features 
#' @return plot a plot object from ggplot2
#' @export
corr_heatmap_fvl_ggplot2 <- function(global_scores,
                                     feature_mat, 
                                     feature_name="Feature", 
                                     select_features=NULL){
  
    # correlate the factors and task
    lf_corrs = cor(global_scores, feature_mat)
    lf_corrs = melt(lf_corrs)
    colnames(lf_corrs) = c('Factor', 'Feature', 'Corr')
    
    # make a heatmap of the correlations
    title = sprintf('Factors versus %s', feature_name)
    plot = ggplot(lf_corrs, aes(Feature, Factor, fill=Corr)) + 
        geom_tile() + 
        ggtitle(title) + 
        theme(plot.title = element_text(size=12)) + 
        scale_fill_distiller(palette = 'RdYlBu', limits=c(-1, 1)) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        scale_y_continuous(breaks=seq(0,nrow(lf_corrs),1))
    
    return(plot)
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
    
    # Annotation databases used for biological enrichment
    path.database <- "data/bio_annotations/c2.cp.reactome.v6.2.symbols.gmt" #REACTOME

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

=======
library("ggplot2")
library("reshape2")
library("ComplexHeatmap")
library("circlize")

#library("dplyr")
#library('biomaRt')

#' Plotting a heatmap of global_loadings versus features
#'
#' @description Plots a heatmap of MCIA global_loadings versus factors
#' @param global_loadings the global_loadings matrix after running MCIA
#' @param data_blocks the list of matrices used to calculate the MCIA factors
#' @param omic_name the name of the omic that should be plot, should be in data_blocks
#' @param select_features a vector of numbers to filter features 
#' @return the ggplot2 object
#' @export
global_loadings_heatmap_ComplexHeatmap <- function(global_loadings,
                                                   data_blocks,
                                                   omic_name, 
                                                   select_features=NULL){
    
        # get the index of the omics we need
        omics_index = which(names(data_blocks) == omic_name)
        
        # find at what index we should start in the global loadings
        if (omics_index == 1){
            global_start = 1 
        } else {
            i = 1
            global_start = 0
            while (i != omics_index){
                global_start = global_start + ncol(data_blocks[[i]])
                i = i + 1
            }
        }
        
        # global_end is global_start plus the number or rows in the 
        # data_blocks of omic_name
        global_end = global_start + ncol(data_blocks[[omics_index]]) - 1
        
        # extract data for current omic
        # transpose so rows are the latent factors and columns are the features
        omic_data = t(global_loadings[seq(global_start, global_end),])
        
        # rename rows and columns for plotting
        rownames(omic_data) = paste0("LF", seq(1, nrow(omic_data)))
        colnames(omic_data) = colnames(data_blocks[[omic_name]])
        
        # make a heatmap of the correlations
        #color_func = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
        coltitle = sprintf(sprintf("%s Features", omic_name))
        
        # filter features basd on select_features 
        if (!is.null(select_features)){
            omic_data = omic_data[, select_features]
        }
        
        # plot the heatmap
        p = Heatmap(omic_data, 
                name = "GL Score", 
                column_title = coltitle,
                row_title = "Latent Factors",
                row_names_gp = gpar(fontsize = 7),
                show_column_names = T,
                show_row_names = T,
                row_names_side = "right"
        )    
      return(p)
}



#' Plotting a heatmap of global_loadings versus features
#'
#' @description Plots a heatmap of MCIA global_loadings versus factors
#' @param global_loadings the global_loadings matrix after running MCIA
#' @param data_blocks the list of matrices used to calculate the MCIA factors
#' @param omic_name the name of the omic that should be plot, should be in data_blocks
#' @param select_features a vector of numbers to filter features 
#' @return the ggplot2 object
#' @export
global_loadings_heatmap_ComplexHeatmap <- function(global_loadings,
                                                   data_blocks,
                                                   omic_name, 
                                                   select_features=NULL){
    
        # get the index of the omics we need
        omics_index = which(names(data_blocks) == omic_name)
        
        # find at what index we should start in the global loadings
        if (omics_index == 1){
            global_start = 1 
        } else {
            i = 1
            global_start = 0
            while (i != omics_index){
                global_start = global_start + ncol(data_blocks[[i]])
                i = i + 1
            }
        }
        
        # global_end is global_start plus the number or rows in the 
        # data_blocks of omic_name
        global_end = global_start + ncol(data_blocks[[omics_index]]) - 1
        
        # extract data for current omic
        # transpose so rows are the latent factors and columns are the features
        omic_data = t(global_loadings[seq(global_start, global_end),])
        
        # rename rows and columns for plotting
        rownames(omic_data) = paste0("LF", seq(1, nrow(omic_data)))
        colnames(omic_data) = colnames(data_blocks[[omic_name]])
        
        # make a heatmap of the correlations
        #color_func = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
        coltitle = sprintf(sprintf("%s Features", omic_name))
        
        # filter features basd on select_features 
        if (!is.null(select_features)){
            omic_data = omic_data[, select_features]
        }
        
        # plot the heatmap
        p = Heatmap(omic_data, 
                name = "GL Score", 
                column_title = coltitle,
                row_title = "Latent Factors",
                row_names_gp = gpar(fontsize = 7),
                show_column_names = T,
                show_row_names = T,
                row_names_side = "right"
        )    
      return(p)
}



#' Plot a correlation heatmap of global scores versus features using ComplexHeatmap
#'
#' @description Plots a heatmap of MCIA factos versus other vectors such as
#' the original feature vectors (membership analysis) or versus tasks to 
#' investigate the predictive power of these factors to some other vector of data.
#' @param global_scores the MCIA factor object calculated from block_matrix
#' @param feature_mat the input matrix used to calculate the MCIA factors
#' @param feature_name the name of the input
#' @param select_features a vector of numbers to filter features 
#' @return the ggplot2 object
#' @export
corr_heatmap_fvl_ComplexHeatmap <- function(global_scores,
                                            feature_mat,
                                            feature_name="Feature",
                                            select_features=NULL){
        
        # correlate the factors and task
        lf_corrs = cor(global_scores, feature_mat)
        rownames(lf_corrs) = paste("LF", seq(1, nrow(lf_corrs)))
        
        # filter features basd on select_features 
        if (!is.null(select_features)){
            lf_corrs = lf_corrs[, select_features]
        }
        
        # make a heatmap of the correlations
        color_func = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
        title = sprintf('MCIA latent factors versus %s', feature_name)
        p = Heatmap(lf_corrs, 
                name = "Pearson's R", 
                column_title = title,
                row_title = "Latent Factors",
                row_names_gp = gpar(fontsize = 7),
                col = color_func,
                show_column_names = T,
                show_row_names = T,
                row_names_side = "right"
        )    
      return(p)
}


#' Plot a correlation heatmap of global scores versus features using ggplot2
#'
#' @description Plots a heatmap of MCIA factos versus other vectors such as
#' the original feature vectors (membership analysis) or versus tasks to 
#' investigate the predictive power of these factors to some other vector of data.
#' @param global_scores the MCIA factor object calculated from block_matrix
#' @param feature_mat the input matrix used to calculate the MCIA factors
#' @param feature_name the name of the input
#' @param select_features a vector of numbers to filter features 
#' @return the ggplot2 object
#' @export
corr_heatmap_fvl_ggplot2 <- function(global_scores,
                                     feature_mat, 
                                     feature_name="Feature", 
                                     select_features=NULL){
  
    # correlate the factors and task
    lf_corrs = cor(global_scores, feature_mat)
    lf_corrs = melt(lf_corrs)
    colnames(lf_corrs) = c('Factor', 'Feature', 'Corr')
    
    # make a heatmap of the correlations
    title = sprintf('MCIA latent factors versus %s', feature_name)
    p = ggplot(lf_corrs, aes(Feature, Factor, fill=Corr)) + 
        geom_tile() + 
        ggtitle(title) + 
        theme(plot.title = element_text(size=12)) + 
        scale_fill_distiller(palette = 'RdYlBu', limits=c(-1, 1)) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        scale_y_continuous(breaks=seq(0,nrow(lf_corrs),1))
    
    return(p)
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
    
    # Annotation databases used for biological enrichment
    path.database <- "data/bio_annotations/c2.cp.reactome.v6.2.symbols.gmt" #REACTOME

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

>>>>>>> fbe14f63785642c9da43fd4de5c1a18c40d17ad7
