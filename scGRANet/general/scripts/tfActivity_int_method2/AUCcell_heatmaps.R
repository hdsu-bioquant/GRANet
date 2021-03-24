options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

p_exprs              <- as.character(args[1])
p_regulonAUC         <- as.character(args[2])
p_regulonAUCbin      <- as.character(args[3])
p_regulonAUC_devCor  <- as.character(args[4])
p_regulonAUCpySCENIC <- as.character(args[5])
p_regulonAUCHeat     <- as.character(args[6])
p_regulonAUCbinHeat  <- as.character(args[7])
p_regulonAUCCorHeat  <- as.character(args[8])
p_regulonAUCpyHeat   <- as.character(args[9])

#------------------------------------------------------------------------------#
p_exprs              <- 'results/scrna/seurat/rna_seurat_int_annot_transfer.RDS'
p_regulonAUC         <- "results/integrated/TF_activity_method2/regulonAUC.RDS"
# p_regulonAUCbin      <- "results/integrated/TF_activity/regulonAUCbin.RDS"
p_regulonAUC_devCor  <- "results/integrated/TF_activity_method2/regulonAUC_devCorrected.RDS"
# p_regulonAUCpySCENIC <- "results/scrna/SCENIC/auc_mtx.csv"
p_regulonAUCHeat     <- "figures/tfActivity_integrated/method2_regulonAUC.pdf"
# p_regulonAUCbinHeat  <- "figures/tfActivity_integrated/regulonAUCbin.pdf"
p_regulonAUCCorHeat  <- "figures/tfActivity_integrated/method2_regulonAUC_devCorrected.pdf"
# p_regulonAUCpyHeat   <- "figures/tfActivity_integrated/regulonAUC_pySCENIC.pdf"
#------------------------------------------------------------------------------#

library(tidyverse)
library(Seurat)
library(ComplexHeatmap)
library(viridis)
library(Polychrome)

##-----------------------------------------------------------------------------#
##                        Regulon heatmap annotation                           #
##-----------------------------------------------------------------------------#
exprs         <- readRDS(p_exprs)
regulonAUC    <- readRDS(p_regulonAUC)
#regulonAUCbin <- readRDS(p_regulonAUCbin)
regulonAUC_devCorrected <- readRDS(p_regulonAUC_devCor)


Heatmap_Regulon <- function(regulon, seurat, normRows, title){
  
  type.colVector <- seurat@meta.data[,"predicted.id", drop=FALSE] %>% 
    dplyr::distinct() %>% 
    dplyr::rename(Celltype =  predicted.id) %>% 
    dplyr::arrange(Celltype) %>% 
    dplyr::mutate(color = alphabet.colors(n())) %>% 
    deframe()
  
  # Build Heatmap annotation
  heat.anno <- HeatmapAnnotation(df  = seurat@meta.data[colnames(regulon),"predicted.id", drop=FALSE],
                                 col = list(predicted.id = type.colVector),
                                 show_annotation_name = TRUE, na_col = "white")
  if (normRows) {
    x <- regulon/rowMaxs(regulon)
  } else {
    x <- regulon
  }
  
  Heatmap(x,
          col = viridis(100),
          name = title,
          top_annotation = heat.anno,
          show_column_names = FALSE,
          show_row_names = FALSE,
          show_column_dend = FALSE,
          show_row_dend = FALSE,
          cluster_rows = TRUE)
}


#------------------------------------------------------------------------------#
#                                Heatmap regulon AUC                           #
#------------------------------------------------------------------------------#
dim(regulonAUC)
p <- Heatmap_Regulon(regulon=regulonAUC, seurat=exprs, normRows=TRUE,
                     title = "Regulon AUC\nmethod 2")
pdf(p_regulonAUCHeat, width = 10, height = 6)
p
dev.off()


p <- Heatmap_Regulon(regulon=regulonAUC_devCorrected, seurat=exprs, normRows=TRUE,
                     title = "Regulon AUC\nmethod 2\ncorrected\nMotif deviation")
pdf(p_regulonAUCCorHeat, width = 10, height = 6)
p
dev.off()


regulonStats <- function(regulon){
  message(paste0("No. Cell: ", ncol(regulon)))
  message(paste0("Regulons: ", nrow(regulon)))
  message(paste0("Positive Regulated: ", length(grep("\\(\\+\\)", rownames(regulon)))))
  message(paste0("Negative Regulated: ", length(grep("\\(\\-\\)", rownames(regulon)))))
}

regulonStats(regulonAUC)
regulonStats(regulonAUC_devCorrected)


