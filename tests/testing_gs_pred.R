## Testing gs prediction 

install.packages("C:/Users/maxim/Documents/NIPALS_MCIA",
                 repos = NULL,
                 type = "source")
library(nipalsMCIA)

data(NCI60)
# performing MCIA
mcia_results <- nipals_multiblock(data_blocks, preprocMethod='colprofile',num_PCs = 10, plots='all', tol=1e-12)

# Predicting "new" global scores with same data
df <- data_blocks
gs_new <- predict_gs(mcia_results,df)

# Should be close to machine precision:
testdiff <- norm(gs_new - mcia_results$global_scores,type="F")
print(testdiff)

# Bonus - heatmap plot of block score weights
MCIA_plots(mcia_results,'block_weights_heatmap')

# Bonus - coloring different clusters
par(mfrow=c(1,1))
CNS = 1:6; LEU = 7:12; ME = 13:21;
clus_list <- list(CNS, LEU, ME)
MCIA_plots(test,'projection',orders = c(1,2),clusters = clus_list, legend_loc = "topright")


