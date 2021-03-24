options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

p_archrproj          <- dirname(as.character(args[1]))
p_exprs              <- as.character(args[2])
p_regulonAUC         <- as.character(args[3])
p_regulonAUC_devCor  <- as.character(args[4])

#------------------------------------------------------------------------------#
# p_archrproj          <- "results/scatac/archr/Save-projcovid6/"
# p_exprs              <- 'results/scrna/seurat/rna_seurat_int_annot_transfer.RDS'
# p_regulonAUC         <- "results/integrated/TF_activity/regulonAUC.RDS"
# #p_regulonAUCbin      <- "results/integrated/TF_activity/regulonAUCbin.RDS"
# p_regulonAUC_devCor  <- "results/integrated/TF_activity/regulonAUC_devCorrected.RDS"
#------------------------------------------------------------------------------#

library(tidyverse)
library(ArchR)
library(Seurat)

#------------------------------------------------------------------------------#
#                      Get motif deviation from ChromVar                       #
#------------------------------------------------------------------------------#
exprs         <- readRDS(p_exprs)
regulonAUC    <- readRDS(p_regulonAUC)
#regulonAUCbin <- readRDS(p_regulonAUCbin)
archrproj <- loadArchRProject(path = p_archrproj)

# Extract motif deviations
motifDeviation <- ArchR::getMatrixFromProject(archrproj, useMatrix = "MotifMatrix")
motifDeviation <- as.matrix(motifDeviation@assays@data$z)

# Mean deviation by cell type
motifDeviation_grouped <- as.tibble(t(motifDeviation)) %>% 
  mutate(CellID = colnames(motifDeviation)) %>% 
  mutate(cellType = archrproj$predictedGroup[match(CellID, archrproj$cellNames)]) %>% 
  dplyr::select(-CellID) %>% 
  group_by(cellType) %>% 
  summarise_all(mean) %>% 
  column_to_rownames("cellType") %>% 
  t()
  
#------------------------------------------------------------------------------#
#                  Motif deviation matched to regulon AUC                      #
#------------------------------------------------------------------------------#

# Match deviation by cell type with scRNAseq data
RegulonMotif_map <- tibble(RegulonID = rownames(regulonAUC)) %>% 
  mutate(TF = sub(" .*", "", RegulonID)) %>% 
  left_join(tibble(TF      = sub("_.*", "", rownames(motifDeviation_grouped)),
                   motifID = rownames(motifDeviation_grouped)),
            by = "TF")

exprs_fil <- exprs[,exprs$predicted.id %in% colnames(motifDeviation_grouped)]
dim(exprs_fil)
motifDeviationAUC <- motifDeviation_grouped[match(RegulonMotif_map$motifID, rownames(motifDeviation_grouped)),
                                            match(exprs_fil$predicted.id, colnames(motifDeviation_grouped))]
colnames(motifDeviationAUC) <- colnames(exprs_fil)

dim(motifDeviationAUC)
min(motifDeviationAUC)
motifDeviationAUC[motifDeviationAUC < 0] <- 0

#------------------------------------------------------------------------------#
#              Correct regulons with motif deviations by cell type             #
#------------------------------------------------------------------------------#

regulonAUC_devCorrected <- regulonAUC[,colnames(motifDeviationAUC)] * motifDeviationAUC
saveRDS(regulonAUC_devCorrected, file = p_regulonAUC_devCor)


# grep("CTCF", rownames(regulonAUC_devCorrected), value = TRUE)
# 
# UMAP_1
# p1 <- as.data.frame(exprs_fil@reductions$umap@cell.embeddings) %>% 
#   rownames_to_column("CellID") %>% 
#   mutate(CTCF_posReg = regulonAUC_devCorrected["CTCF (+)", match(CellID, colnames(regulonAUC_devCorrected))]) %>% 
#   mutate(CTCF_posReg = CTCF_posReg/max(CTCF_posReg)) %>% 
#   ggplot(aes(x=UMAP_1, y=UMAP_2, color=CTCF_posReg)) +
#   geom_point() +
#   scale_color_viridis() +
#   cowplot::theme_cowplot()
# 
# p2 <- as.data.frame(exprs_fil@reductions$umap@cell.embeddings) %>% 
#   rownames_to_column("CellID") %>% 
#   mutate(CellType = exprs_fil$predicted.id) %>% 
#   ggplot(aes(x=UMAP_1, y=UMAP_2, color=CellType)) +
#   geom_point() +
#   cowplot::theme_cowplot()
# 
# 
# library(patchwork)
# p1+p2






# colnames(motifDeviationAUC)
# 
# #rowMaxs(motifDeviationAUC/rowMaxs(motifDeviationAUC))
# 
# 
# 
# 
# 
# seATAC <- readRDS("results/scatac/archr/Save-projcovid6/Peak2GeneLinks/seATAC-Group-KNN.rds")
# seRNA  <- readRDS("results/scatac/archr/Save-projcovid6/Peak2GeneLinks/seRNA-Group-KNN.rds")
# length(seATAC@metadata$KNNList)
# length(seRNA@metadata$KNNList)
# 
# length(unique(unlist(seRNA@metadata$KNNList)))
# length(unique(unlist(seATAC@metadata$KNNList)))
# 
# 
# 
# z <- tibble(groupedID  = unique(unlist(seRNA@metadata$KNNList))) %>% 
#   mutate(CommonID = sub(".*#", "", groupedID)) %>% 
#   left_join(tibble(OriginalID = colnames(exprs),
#                    CommonID   = sub("_.*", "", colnames(exprs))), 
#             by="CommonID")
# 
# 
# 
# dim(exprs)
# 
# x <- ArchR::getMatrixFromProject(archrproj, useMatrix = "GeneIntegrationMatrix")
# x
# 
# table(colnames(x) %in% unique(unlist(seATAC@metadata$KNNList)))
# 
# y <- colnames(x)[!colnames(x) %in% unique(unlist(seATAC@metadata$KNNList))]
# 
# table(y %in% z$CommonID)
# 
# 
#   filter(OriginalID %in% colnames(exprs))
# 
# 
# 
#   
#   
#   
# 
# head(unique(unlist(seRNA@metadata$KNNList)))
# head(colnames(exprs))
# 
# table(colnames(exprs) %in% unique(unlist(seRNA@metadata$KNNList)))
# 
#       dim(exprs)
# 
# as.list(x@metadata$KNNList)
# 
# 
# 
# ArchR::getAvailableMatrices(archrproj)
# ArchR::get(archrproj)
# getPeak2GeneLinks
# 
# addPeak2GeneLinks
# 
# archrproj
# 
# 
# 
# 
# 
# 
# 
# # getPositions()
# # peakAnnoEnrichment
# # x <- ArchR::getMatrixFromProject(archrproj, useMatrix = "MotifMatrix")
# # x@assays@data$deviations
# # rownames(x)
# # colnames(x)
# # 
# # ArchR::getAvailableMatrices(archrproj)
# # ArchR::getPeakAnnotation(archrproj)
# # ArchR::getPeakSet(archrproj)
# # 
# # getMatches(archrproj)
# # # ChromVar motif deviations
# # dfVarDev <- getVarDeviations(archrproj, name = "MotifMatrix", plot = FALSE)
# # dim(dfVarDev)
# # #length(unique(dfVarDev$idx))
# # 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# exprs         <- readRDS(p_exprs)
# regulonAUC    <- readRDS(p_regulonAUC)
# regulonAUCbin <- readRDS(p_regulonAUCbin)
# 
# 
# pal2 <- alphabet.colors(26)
# ranpoints(pal2, 14)
# 
# 
# #exprs@meta.data[colnames(regulonAUC),"predicted.id", drop=FALSE] %>% 
# type.colVector <- exprs@meta.data[,"predicted.id", drop=FALSE] %>% 
#   dplyr::distinct() %>% 
#   dplyr::rename(Celltype =  predicted.id) %>% 
#   dplyr::arrange(Celltype) %>% 
#   dplyr::mutate(color = alphabet.colors(n())) %>% 
#   deframe()
# 
# 
# #------------------------------------------------------------------------------#
# #                          Heatmap binarized regulon AUC                       #
# #------------------------------------------------------------------------------#
# # Build Heatmap annotation
# heat.anno <- HeatmapAnnotation(df  = exprs@meta.data[colnames(regulonAUC),"predicted.id", drop=FALSE],
#                                col = list(predicted.id = type.colVector),
#                                show_annotation_name = TRUE, na_col = "white")
# 
# p <- Heatmap(regulonAUC,
#              col = viridis(100),
#              name = paste0("Regulon AUC"),
#              #clustering_distance_columns = 'pearson',
#              top_annotation = heat.anno,
#              show_column_names = FALSE,
#              show_row_names = FALSE,
#              show_column_dend = FALSE,
#              show_row_dend = FALSE,
#              cluster_rows = TRUE)
# #p
# 
# pdf(p_regulonAUCHeat, width = 10, height = 6)
# p
# dev.off()
# 
# #------------------------------------------------------------------------------#
# #                          Heatmap binarized regulon AUC                       #
# #------------------------------------------------------------------------------#
# 
# # Build Heatmap annotation
# heat.anno <- HeatmapAnnotation(df  = exprs@meta.data[colnames(regulonAUCbin),"predicted.id", drop=FALSE],
#                                col = list(predicted.id = type.colVector),
#                                show_annotation_name = TRUE, na_col = "white")
# 
# p <- Heatmap(regulonAUCbin,
#              col = viridis(100),
#              name = paste0("Regulon AUC"),
#              #clustering_distance_columns = 'pearson',
#              top_annotation = heat.anno,
#              show_column_names = FALSE,
#              show_row_names = FALSE,
#              show_column_dend = FALSE,
#              show_row_dend = FALSE,
#              cluster_rows = TRUE)
# #p
# 
# pdf(p_regulonAUCbinHeat, width = 10, height = 6)
# p
# dev.off()
# 
# 
