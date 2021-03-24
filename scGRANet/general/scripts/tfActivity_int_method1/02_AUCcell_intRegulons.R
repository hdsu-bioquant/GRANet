options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

p_exprs          <- as.character(args[1])
p_regulons       <- as.character(args[2])
p_aucellRankings <- as.character(args[3])
p_regulonAUC     <- as.character(args[4])
#p_regulonAUCbin  <- as.character(args[5])
nCores           <- as.numeric(args[5])

#------------------------------------------------------------------------------#
# p_exprs          <- 'results/scrna/seurat/rna_seurat_int_annot_transfer.RDS'
# p_regulons       <- "results/integrated/TF_activity/tfRegulons_asDF.RDS"
# p_aucellRankings <- "results/integrated/TF_activity/aucellRankings.RDS"
# p_regulonAUC     <- "results/integrated/TF_activity/regulonAUC.RDS"
# p_regulonAUCbin  <- "results/integrated/TF_activity/regulonAUCbin.RDS"
# nCores           <- 20
#------------------------------------------------------------------------------#

library(tidyverse)
library(Seurat)
library(AUCell)
#library(ArchR)


#------------------------------------------------------------------------------#
#               Read Expression matrix and make rankings.                      #
#------------------------------------------------------------------------------#
# Create rankings
exprs <- readRDS(p_exprs)
aucellRankings <- AUCell_buildRankings(exprs@assays$RNA@counts, nCores=nCores, 
                                       plotStats=FALSE, verbose=TRUE)

# Save expression rankings
saveRDS(aucellRankings, file=p_aucellRankings)

#------------------------------------------------------------------------------#
#                               Regulons as list                               #
#------------------------------------------------------------------------------#
regulons_df <- readRDS(p_regulons)

# Function to split df by groups into a named list
named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = " / ")))
  
  grouped %>% 
    dplyr::group_split() %>% 
    rlang::set_names(names)
}

# make list of regulons
regulons_list <- regulons_df %>% 
  dplyr::mutate(RegulonID = paste0(TF, " (",
                            if_else(regulation == 1, '+', "-"),
                            ")")) %>% 
  named_group_split(RegulonID) %>% 
  # add TF to regulon
  purrr::map(function(x) c(unique(x$TF), x$target))



#------------------------------------------------------------------------------#
#               Calculate AUC for each regulon in each cell                    #
#------------------------------------------------------------------------------#
# Calculate AUC
regulonAUC <- AUCell_calcAUC(regulons_list, aucellRankings, 
                             aucMaxRank = aucellRankings@nGenesDetected["1%"], 
                             nCores     = nCores)

# Order the modules by similarity, for easier exploration in the upcoming steps & save
regulonOrder <- orderAUC(regulonAUC) # added to AUCell 1.5.1
regulonAUC <- regulonAUC[regulonOrder,]
regulonAUC

saveRDS(regulonAUC@assays@data$AUC, file=p_regulonAUC)

# #------------------------------------------------------------------------------#
# #                           Binarize regulon activity                          #
# #------------------------------------------------------------------------------#
# cells_AUCellThresholds <- AUCell_exploreThresholds(regulonAUC, 
#                                                    #smallestPopPercent=getSettings(scenicOptions,"aucell/smallestPopPercent"),
#                                                    assignCells=TRUE, plotHist=FALSE, 
#                                                    verbose=FALSE, nCores=1)
# thresholds <- getThresholdSelected(cells_AUCellThresholds)
# 
# 
# # Assign cells
# regulonsCells <- setNames(lapply(names(thresholds), 
#                                  function(x) {
#                                    trh <- thresholds[x]
#                                    names(which(getAUC(regulonAUC)[x,]>trh))
#                                  }),names(thresholds))
# ### Convert to matrix (regulons with zero assigned cells are lost)
# regulonActivity <- reshape2::melt(regulonsCells)
# binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
# class(binaryRegulonActivity) <- "matrix"
# 
# saveRDS(binaryRegulonActivity, file=p_regulonAUCbin)
# 
