options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

p_exprs      <- as.character(args[1])
p_regulons   <- as.character(args[2])
p_regulonAUC <- as.character(args[3])
nCores       <- as.numeric(args[4])
annotCol     <- as.character(args[5])
#------------------------------------------------------------------------------#
# p_exprs          <- 'results/scrna/seurat/rna_seurat_int_annot_transfer.RDS'
# p_regulons       <- "results/integrated/TF_activity_method2/tfRegulons_asDF_byCell.RDS"
# p_regulonAUC     <- "results/integrated/TF_activity_method2/regulonAUC.RDS"
# p_regulonAUCbin  <- "results/integrated/TF_activity_method2/regulonAUCbin.RDS"
# nCores           <- 20
# p_exprs <-  "~/projects/mouse_atlas_TF_activity/data/scrna/Han_atlas/MCA_AdultLung_Seurat.RDS"
# annotCol     <- "Annotation"

# p_exprs          <- '/home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scrna/Han_atlas/MCA_AdultLung_Seurat.RDS'
# p_regulons       <- "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_lung/results/integrated/TF_activity_method2_5000/tfRegulons_asDF.RDS"
# p_regulonAUC     <- "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_lung/results/integrated/TF_activity_method2_5000/regulonAUC.RDS"
# nCores           <- 20
# annotCol     <- "Annotation"

#------------------------------------------------------------------------------#

library(tidyverse)
library(Seurat)
library(AUCell)
#library(ArchR)

#head(Seurat::FetchData(exprs, vars=annotCol)[,1])
#------------------------------------------------------------------------------#
#               Read Expression matrix and make rankings.                      #
#------------------------------------------------------------------------------#
# Create rankings
exprs <- readRDS(p_exprs)
regulons_df_Cell <- readRDS(p_regulons)

# Filter cell for which cell type specific regulon is not available
#exprs <- exprs[,exprs$predicted.id %in% names(regulons_df_Cell)]
exprs <- exprs[,Seurat::FetchData(exprs, vars=annotCol)[,1] %in% names(regulons_df_Cell)]


aucellRankings <- AUCell_buildRankings(exprs@assays$RNA@counts, nCores=nCores, 
                                       plotStats=FALSE, verbose=TRUE)
dim(aucellRankings)

#------------------------------------------------------------------------------#
#                               Regulons as list                               #
#------------------------------------------------------------------------------#
# Function to split df by groups into a named list
named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = " / ")))
  
  grouped %>% 
    dplyr::group_split() %>% 
    rlang::set_names(names)
}

# make list of regulons
regulons_list_Cell <- lapply(regulons_df_Cell, function(regulons_df){
  regulons_df %>% 
    dplyr::mutate(RegulonID = paste0(TF, " (",
                                     if_else(regulation == 1, '+', "-"),
                                     ")")) %>% 
    named_group_split(RegulonID) %>% 
    # add TF to regulon
    purrr::map(function(x) c(unique(x$TF), x$target))
  
})

#------------------------------------------------------------------------------#
#               Calculate AUC for each regulon in each cell                    #
#------------------------------------------------------------------------------#
# Calculate AUC
regulonAUC_Cell <- lapply(setNames(names(regulons_list_Cell), names(regulons_list_Cell)), function(CellType){
  #CellIDs <- colnames(exprs)[exprs$predicted.id %in% CellType]
  CellIDs <- colnames(exprs)[Seurat::FetchData(exprs, vars=annotCol)[,1] %in% CellType]
  
  AUCell_calcAUC(regulons_list_Cell[[CellType]], aucellRankings[,CellIDs], 
                 aucMaxRank = aucellRankings@nGenesDetected["1%"], 
                 nCores     = nCores)
})

regulonAUC_Cell_df <- lapply(regulonAUC_Cell, function(regulonAUC_byCell){
  #dim(regulonAUC_byCell)
  as.data.frame(regulonAUC_byCell@assays@data$AUC) %>% 
    rownames_to_column("Regulon") 
})


regulonAUC <- regulonAUC_Cell_df %>% 
  purrr::reduce(full_join, by = "Regulon") %>% 
  column_to_rownames("Regulon") 

dim(regulonAUC)
sum(is.na(regulonAUC))
regulonAUC[is.na(regulonAUC)] <- 0
regulonAUC <- as.matrix(regulonAUC)
regulonAUC <- regulonAUC[rowSums(regulonAUC) > 0,]


saveRDS(regulonAUC, file=p_regulonAUC)



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
