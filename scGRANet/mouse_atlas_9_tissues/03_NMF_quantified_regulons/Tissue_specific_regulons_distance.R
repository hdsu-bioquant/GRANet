params <- list(p_regulonByType ="/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/tfRegulons_asDF.RDS",
               p_regulonDist ="/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/tfRegulons_dist.RDS",
               p_nmfAUCpySCENIC = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/scrna/SCENIC/NMF_whole/NMF_regulonAUC.RDS",
               p_nmfregulonAUC  = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_regulonAUC.RDS",
               K=11)

regulonByType <- readRDS(params$p_regulonByType)

regulonByType <- lapply(names(regulonByType), function(id){
  regulonByType[[id]] %>% 
    mutate(tissue = id)
})


regulonByTypedf <- bind_rows(regulonByType)

regulonBinMatrix <- regulonByTypedf %>% 
  ungroup() %>% 
  mutate(RegulonID = paste0(tissue, "_", TF)) %>% 
  select(RegulonID, target) %>% 
  mutate(p = 1) %>% 
  pivot_wider(names_from = target, values_from = p) %>% 
  column_to_rownames("RegulonID") %>% 
  as.matrix()
dim(regulonBinMatrix)
regulonBinMatrix[is.na(regulonBinMatrix)] <- 0
regulonBinMatrix[1:5, 1:5]

#install.packages("parallelDist")
library(parallelDist)
#regulonBinDist <- parallelDist::parDist(regulonBinMatrix, method="binary", threads=50)
#as.matrix(regulonBinDist)
#saveRDS(regulonBinDist, params$p_regulonDist)
regulonBinDist <- readRDS(params$p_regulonDist)

regulonBinDistMat <- as.matrix(regulonBinDist)
regulonBinDistMat[1:5,1:5]


colnames(regulonBinDistMat)

#------------------------------------------------------------------------------#
#                     Signature specific regulons                              #
#------------------------------------------------------------------------------#

# Search Brain signature
# Signature 6 for method 2
# Signature 1 for pySCENIC

# Search Thymus signature
# Signature 10 for method 2
# Signature 3 for pySCENIC

# Search kidney signature
# Signature 9 for method 2
# Signature 10 for pySCENIC

# Search liver signature
# Signature 2 for method 2
# Signature 5 for pySCENIC

SSFmintRegul <- SignatureSpecificFeatures(nmfregulonAUC, k=params$K, return_all_features=TRUE)
SSFmpySCENIC <- SignatureSpecificFeatures(nmfAUCpySCENIC, k=params$K, return_all_features=TRUE)

SSFmintRegul <- SSFmintRegul[,c(6,10,9,2)]
SSFmpySCENIC <- SSFmpySCENIC[,c(1,3,10,5)]
colnames(SSFmintRegul) <- colnames(SSFmpySCENIC) <- c("Brain", "Thymus", "Kidney", "Liver")

activeRegulons <- as.data.frame(SSFmintRegul) %>% 
  rownames_to_column("RegulonID") %>% 
  dplyr::filter(!grepl(" \\(-.*", RegulonID)) %>% 
  mutate(RegulonID = sub(" \\(.*", "", RegulonID)) %>% 
  pivot_longer(names_to = "tissue", cols = -RegulonID) %>% 
  dplyr::filter(value == 1) %>% 
  mutate(TF = RegulonID) %>% 
  mutate(RegulonID = paste0("Adult", tissue, "_", TF)) 

#------------------------------------------------------------------------------#
#                     Distance active regulons.                                #
#------------------------------------------------------------------------------#
actfil <- activeRegulons %>% 
  group_by(TF) %>% 
  dplyr::filter(n()<=2)


activeRegulonsDist <- regulonBinDistMat[rownames(regulonBinDistMat) %in% actfil$RegulonID,
                                        colnames(regulonBinDistMat) %in% actfil$RegulonID]

activeRegulonsDist <- 1-activeRegulonsDist
dim(activeRegulonsDist)
diag(activeRegulonsDist) <- 1
heatAnnot <- HeatmapAnnotation(df=as.data.frame(actfil[match(colnames(activeRegulonsDist), actfil$RegulonID),"tissue"]))
Heatmap(activeRegulonsDist,
        col = circlize::colorRamp2(seq(from=0, to=quantile(activeRegulonsDist, .99), length.out=100),
                                   inferno(100)),
        top_annotation = heatAnnot,
        name = "Regulon distance",
        use_raster = TRUE,
        raster_device = "png",
        show_column_names = FALSE,
        show_row_names = FALSE,
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        cluster_rows = TRUE)
# 


#------------------------------------------------------------------------------#
#                     Specific regulon for each tissue                         #
#------------------------------------------------------------------------------#

# Search Brain signature          # Signature 6 for method 2
# Search Thymus signature         # Signature 10 for method 2
# Search kidney signature         # Signature 9 for method 2
# Search liver signature          # Signature 2 for method 2
# Search bone marrow signature    # Signature 8 for method 2
# Search Lung signature           # Signature 1 for method 2
# Search SmallIntestine signature # Signature 3 for method 2
# Search Spleen signature         # Signature 7 for method 2
# Search Testis signature         # Signature 4 for method 2

SSFmintRegul <- SignatureSpecificFeatures(nmfregulonAUC, k=params$K, return_all_features=TRUE)
SSFmintRegul <- SSFmintRegul[,c(6,10,9,2,8,1,3,7,4)]
colnames(SSFmintRegul) <- c("Brain", "Thymus", "Kidney", "Liver", "BoneMarrow",
                            "Lung", "SmallIntestine", "Spleen", "Testis")

activeRegulons <- as.data.frame(SSFmintRegul) %>% 
  rownames_to_column("RegulonID") %>% 
  dplyr::filter(!grepl(" \\(-.*", RegulonID)) %>% 
  mutate(RegulonID = sub(" \\(.*", "", RegulonID)) %>% 
  pivot_longer(names_to = "Tissue", cols = -RegulonID) %>% 
  dplyr::filter(value == 1) %>% 
  mutate(TF = RegulonID) %>% 
  mutate(RegulonID = paste0("Adult", Tissue, "_", TF)) 

activeRegulons


tissue_colors <- setNames(alphabet.colors(9),  c("BoneMarrow", "Brain", "Kidney", "Liver", "Lung", 
                                "SmallIntestine", "Spleen", "Testis", "Thymus"))


#------------------------------------------------------------------------------#
#                     Distance general regulons.                               #
#------------------------------------------------------------------------------#
heat_list <- list()

actfil <- activeRegulons %>% 
  group_by(TF) %>% 
  dplyr::filter(n()>=6)
actfil

table(actfil$TF)

activeRegulonsDist <- regulonBinDistMat[rownames(regulonBinDistMat) %in% actfil$RegulonID,
                                        colnames(regulonBinDistMat) %in% actfil$RegulonID]

activeRegulonsDist <- 1-activeRegulonsDist
dim(activeRegulonsDist)
diag(activeRegulonsDist) <- 1
heatAnnot <- HeatmapAnnotation(df=as.data.frame(actfil[match(colnames(activeRegulonsDist),
                                                             actfil$RegulonID),c("Tissue", "TF")]),
                               col = list(Tissue = tissue_colors))
heat_list[[1]] <- Heatmap(activeRegulonsDist,
                          col = circlize::colorRamp2(seq(from=0, 
                                                         to=quantile(activeRegulonsDist, .99), 
                                                         length.out=100),
                                                     inferno(100)),
                          top_annotation = heatAnnot,
                          name = "General-Regulons\n(active in 6 or more tissues)\nsimilarity",
                          use_raster = TRUE,
                          raster_device = "png",
                          show_column_names = FALSE,
                          show_row_names = FALSE,
                          show_column_dend = TRUE,
                          show_row_dend = TRUE,
                          cluster_rows = TRUE)
heat_list[[1]]

#------------------------------------------------------------------------------#
#                     Distance specific regulons.                              #
#------------------------------------------------------------------------------#

actfil <- activeRegulons %>% 
  group_by(TF) %>% 
  dplyr::filter(n()==1)
actfil

activeRegulonsDist <- regulonBinDistMat[rownames(regulonBinDistMat) %in% actfil$RegulonID,
                                        colnames(regulonBinDistMat) %in% actfil$RegulonID]

activeRegulonsDist <- 1-activeRegulonsDist
dim(activeRegulonsDist)
#diag(activeRegulonsDist) <- 1
heatAnnot <- HeatmapAnnotation(df=as.data.frame(actfil[match(colnames(activeRegulonsDist),
                                                             actfil$RegulonID),c("Tissue", "TF")]),
                               col = list(Tissue = tissue_colors))
heat_list[[2]] <- Heatmap(activeRegulonsDist,
                          col = circlize::colorRamp2(seq(from=0, 
                                                         to=quantile(activeRegulonsDist, .99), 
                                                         length.out=100),
                                                     inferno(100)),
                          top_annotation = heatAnnot,
                          name = "Specific-Regulons\n(Differentially active in one tissue)\nsimilarity",
                          use_raster = TRUE,
                          raster_device = "png",
                          show_column_names = FALSE,
                          show_row_names = FALSE,
                          show_column_dend = TRUE,
                          show_row_dend = TRUE,
                          cluster_rows = TRUE)
heat_list[[2]]



#------------------------------------------------------------------------------#
#                                 Jun Fos                                      #
#------------------------------------------------------------------------------#

allRegulons <- as.data.frame(SSFmintRegul) %>% 
  rownames_to_column("RegulonID") %>% 
  dplyr::filter(!grepl(" \\(-.*", RegulonID)) %>% 
  mutate(RegulonID = sub(" \\(.*", "", RegulonID)) %>% 
  pivot_longer(names_to = "Tissue", cols = -RegulonID) %>% 
  #dplyr::filter(value == 1) %>% 
  mutate(TF = RegulonID) %>% 
  mutate(RegulonID = paste0("Adult", Tissue, "_", TF)) 

allRegulons


actfil <- allRegulons %>% 
  group_by(TF) %>% 
  dplyr::filter(TF %in% c("Jun", "Fos"))
actfil

activeRegulonsDist <- regulonBinDistMat[rownames(regulonBinDistMat) %in% actfil$RegulonID,
                                        colnames(regulonBinDistMat) %in% actfil$RegulonID]

activeRegulonsDist <- 1-activeRegulonsDist
dim(activeRegulonsDist)
#diag(activeRegulonsDist) <- 1
heatAnnot <- HeatmapAnnotation(df=as.data.frame(actfil[match(colnames(activeRegulonsDist),
                                                             actfil$RegulonID),c("Tissue", "TF")]),
                               col = list(Tissue = tissue_colors))

heat_list[[3]] <- Heatmap(activeRegulonsDist,
                          col = circlize::colorRamp2(seq(from=0, 
                                                         to=quantile(activeRegulonsDist, .99), 
                                                         length.out=100),
                                                     inferno(100)),
                          top_annotation = heatAnnot,
                          name = "Jun/Fos\nRegulon similarity",
                          use_raster = TRUE,
                          raster_device = "png",
                          show_column_names = FALSE,
                          show_row_names = FALSE,
                          show_column_dend = TRUE,
                          show_row_dend = TRUE,
                          cluster_rows = TRUE)
heat_list[[3]]


#------------------------------------------------------------------------------#
#                     Distance specific regulons.                              #
#------------------------------------------------------------------------------#

allRegulons <- as.data.frame(SSFmintRegul) %>% 
  rownames_to_column("RegulonID") %>% 
  dplyr::filter(!grepl(" \\(-.*", RegulonID)) %>% 
  mutate(RegulonID = sub(" \\(.*", "", RegulonID)) %>% 
  pivot_longer(names_to = "Tissue", cols = -RegulonID) %>% 
  #dplyr::filter(value == 1) %>% 
  mutate(TF = RegulonID) %>% 
  mutate(RegulonID = paste0("Adult", Tissue, "_", TF)) 

allRegulons


actfil <- allRegulons %>% 
  group_by(TF) %>% 
  dplyr::filter(sum(value) == 1)
actfil



activeRegulonsDist <- regulonBinDistMat[rownames(regulonBinDistMat) %in% actfil$RegulonID,
                                        colnames(regulonBinDistMat) %in% actfil$RegulonID]

activeRegulonsDist <- 1-activeRegulonsDist
dim(activeRegulonsDist)
#diag(activeRegulonsDist) <- 1
heatAnnot <- HeatmapAnnotation(df=as.data.frame(actfil[match(colnames(activeRegulonsDist),
                                                             actfil$RegulonID),c("Tissue", "TF")]),
                               col = list(Tissue = tissue_colors))
heat_list[[4]] <- Heatmap(activeRegulonsDist,
                          col = circlize::colorRamp2(seq(from=0, 
                                                         to=quantile(activeRegulonsDist, .99), 
                                                         length.out=100),
                                                     inferno(100)),
                          top_annotation = heatAnnot,
                          name = "Specific-Regulons\n(Differentially active in one tissue)\nsimilarity",
                          use_raster = TRUE,
                          raster_device = "png",
                          show_column_names = FALSE,
                          show_row_names = FALSE,
                          show_column_dend = TRUE,
                          show_row_dend = TRUE,
                          cluster_rows = TRUE)
heat_list[[4]]




pdf(file = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_TF_expression_tissue_Regulon_similarity.pdf"
, width=10, height=7)
for (h in heat_list) {
  print(h)
}
dev.off()





