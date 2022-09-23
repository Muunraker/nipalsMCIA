#library('biomaRt')

#' Plotting a heatmap of global factors scores (sample v. factors)
#'
#' @description Plots a heatmap of MCIA global scores 
#' @param global_scores the global_scores matrix after running MCIA
#' @return the ggplot2 object
#' @export
global_scores_heatmap_ComplexHeatmap <- function(global_scores){
    colnames(global_scores) = paste0('F', seq(1, ncol(global_scores)))
    p = ComplexHeatmap::Heatmap(global_scores, 
                                name = "GS Score", 
                                column_title = "Factors",
                                row_title = "Samples",
                                row_names_gp = grid::gpar(fontsize = 7),
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
        p = ComplexHeatmap::Heatmap(omic_data, 
                name = "GL Score", 
                column_title = coltitle,
                row_title = "Latent Factors",
                row_names_gp = grid::gpar(fontsize = 7),
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
        p = ComplexHeatmap::Heatmap(omic_data, 
                name = "GL Score", 
                column_title = coltitle,
                row_title = "Latent Factors",
                row_names_gp = grid::gpar(fontsize = 7),
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
        p = ComplexHeatmap::Heatmap(lf_corrs, 
                name = "Pearson's R", 
                column_title = title,
                row_title = "Latent Factors",
                row_names_gp = grid::gpar(fontsize = 7),
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
#' @description Runs fgsea for the input gene vector
#' @param metagenes Vector of gene scores where the row names are HUGO symbols 
#' @param path.database path to a GMT annotation file
#' @param factor vector of factors which should be analyzed
#' @param pval.thr = p-value threshold (default to 0.05)
#' @return data frame with the most significant p-value number of significant
#' pathways
#' @return the selectivity scores across the given factors
#' @export
gsea_report <- function(metagenes, path.database, factors=NULL, pval.thr=0.05,
                        nproc=4){
    
    # Load annotation database
    pathways <- fgsea::gmtPathways(path.database)
    
    # Number of factors
    if (is.null(factors)){
        factors <- seq(ncol(metagenes))
    }
    
    # Containers to report results
    report_min_pval <- numeric(0)
    report_number <- numeric(0)
    sig_paths <- numeric(0)
    
    # Calculate biological annotation enrichment, for each factor
    for(j in factors){
        
        print(paste0('Running GSEA for Factor', j))
      
        # Extract a vector of scores for GSEA and set rownams
        # to HUGO symbols
        scores <- setNames(as.matrix(metagenes[,j]), rownames(metagenes))
        
        # Compute GSEA
        fgseaRes <- fgsea::fgseaMultilevel(pathways, scores, nPermSimple = 10000,
                                    minSize=15, maxSize=500, nproc=nproc)
    
        
        # Report if at least one pathway is significant and
        # the min-pvalue from all pathways
        if(sum(fgseaRes$padj < pval.thr)!=0){
            
            # Report the minimum p-value
            min_pval = min(fgseaRes$padj)
            report_min_pval <- rbind(report_min_pval, min_pval)
            
            # Keep names of significant pathways
            curr_sig_paths = fgseaRes[fgseaRes$padj < pval.thr, "pathway"]
            sig_paths <- c(sig_paths, curr_sig_paths)
            
            # Report number of unique significant pathways  
            report_number <- rbind(report_number, dplyr::n_distinct(curr_sig_paths))
            
        } else{
            
            # Report the minimum p-value, assigning NA
            report_min_pval <- rbind(report_min_pval, NA)
            
            # Report number of unique significant pathways, assigning 0
            report_number <- rbind(report_number, 0)
        }
    }
    
    # Report selectivity
    # Selectivity is calculated best across ALL factors
    # The formula is give within https://www.nature.com/articles/s41467-020-20430-7
    # and section "Selectivity Score". SS = (Nc + Nf) / 2L
    # Nc is the total number of clinical annotations associated with at
    # least a factor, Nf the total number of factors associated with at
    # least a clinical annotation, and L the total number of
    # associations between clinical annotations and factors. S has a
    # maximum value of 1 when each factor is associated with one and
    # only one clinical/biological annotation, and a minimum of 0 in
    # the opposite case. An optimal method should thus maximize its
    # number of factors associated with clinical/biological annotations
    # without having a too low selectivity.
    Nc <- n_distinct(sig_paths)
    Nf <- length(which(!is.na(min_pval)))
    L <- length(sig_paths)
    SS = (Nc + Nf) / (2* L)
    
    # Setting up the final reporting data structure with list containing
    # two values, the first is a dataframe where the rows are factors and
    # the columns is composed of the following fields:
    # - the minimum p-value of all pathways
    # - the total number of significant pathways
    # Lastly, the second value is the selectivity value which assesses whether
    # the set of factors are capturing very different gene sets/pathways from
    # one another
    out <- list()
    out[[1]] <- data.frame(min_pval=report_min_pval,
                           total_pathways=report_number)
    row.names(out[[1]]) <- paste0('Factor', factors)
    
    out[[2]] <- SS
    names(out) <- c('per-factor-results', 'selectivity')
    return(out)
}