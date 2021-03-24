
params <- list(p_nmfAUCpySCENIC = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/scrna/SCENIC/NMF_whole/NMF_regulonAUC.RDS",
               p_nmfregulonAUC  = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_regulonAUC.RDS",
               K = 11,
               p_cooccurrance = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_TF_tissue_coocurrance.pdf",
               p_exprs = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scrna/Han_atlas/MCA_Adult_9Tissues_Seurat.RDS"
)

#------------------------------------------------------------------------------#
#                               Read Data                                      #
#------------------------------------------------------------------------------#

nmfregulonAUC  <- readRDS(params$p_nmfregulonAUC)
nmfAUCpySCENIC <- readRDS(params$p_nmfAUCpySCENIC)


seuratobj <- readRDS(params$p_exprs)
seuratobj <- NormalizeData(seuratobj, normalization.method = "LogNormalize", scale.factor = 1e6)
seuratobj <- ScaleData(seuratobj, features = rownames(seuratobj))
seuratobj$ident <- sub("Adult", "", seuratobj$tissue)

#------------------------------------------------------------------------------#
#                     Signature specific regulons                              #
#------------------------------------------------------------------------------#

WintRegul <- WMatrix(nmfregulonAUC, k=params$K)
WpySCENIC <- WMatrix(nmfAUCpySCENIC, k=params$K)

regulonID_df <- tibble(Regulon = rownames(WintRegul)) %>% 
  dplyr::filter(!grepl("\\(-\\)", Regulon)) %>% 
  mutate(TF = sub(" \\(.*", "", Regulon)) %>% 
  full_join(tibble(pySCENIC_Regulon = rownames(WpySCENIC),
                   TF = sub("\\(.*", "", rownames(WpySCENIC))),
            bt="TF")
  

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

SSFintRegul <- SignatureSpecificFeatures(nmfregulonAUC, k=params$K)
SSFpySCENIC <- SignatureSpecificFeatures(nmfAUCpySCENIC, k=params$K)

SSFintRegul <- SSFintRegul[c(6,10,9,2)]
SSFpySCENIC <- SSFpySCENIC[c(1,3,10,5)]
names(SSFintRegul) <- names(SSFpySCENIC) <- c("Brain", "Thymus", "Kidney", "Liver")

# Get only name
SSFintRegul <- lapply(SSFintRegul, function(x){
  x <- x[!grepl("\\(-\\)", x)]
  sub(" \\(.*", "", x)
})
SSFpySCENIC <- lapply(SSFpySCENIC, function(x){
  x <- x[!grepl("\\(-\\)", x)]
  sub("\\(.*", "", x)
})



#------------------------------------------------------------------------------#
#                     Expression of each gene by tissue                        #
#------------------------------------------------------------------------------#
#DotPlot(seuratobj, features=SSFintRegul$Brain, group.by="ident") + RotatedAxis()
#DotPlot(seuratobj, features=SSFpySCENIC$Brain, group.by="ident")


uniqueSSF_gg <-  lapply(names(SSFintRegul), function(id){
  ssf <- SSFintRegul[[id]]
  # Get TF identified but our method but not pySCENIC
  x <- regulonID_df[regulonID_df$TF %in% ssf,] %>% 
    dplyr::filter(is.na(pySCENIC_Regulon)) %>% select(TF) %>% deframe()
  
  p <- DotPlot(seuratobj, features=x, group.by="ident", cluster.idents = TRUE) + RotatedAxis() + 
    scale_color_viridis() + 
    xlab(paste0("Integrative regulon specific ", id, " TFs\n(Not found by pySCENIC)")) +
    theme(panel.background = element_rect(fill = "black"))
  return(p)
})
patchwork::wrap_plots(uniqueSSF_gg)

ggsave("/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_TF_expression_tissue_specificIntegrative.pdf", width = 17, height = 7)


x <- regulonID_df[regulonID_df$TF %in% SSFintRegul$Kidney,] %>% 
  dplyr::filter(is.na(pySCENIC_Regulon)) %>% select(TF) %>% deframe()
DotPlot(seuratobj[,seuratobj$ident == "Kidney"], features=x, group.by="Annotation") + 
  RotatedAxis() + scale_color_viridis() + 
  xlab(paste0("Integrative regulon specific Kidney TFs\n(Not found by pySCENIC)")) +
  theme(panel.background = element_rect(fill = "black"))
ggsave("/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_TF_expression_Kidney_specificIntegrative.pdf", width = 9, height = 6)


x <- regulonID_df[regulonID_df$TF %in% SSFintRegul$Thymus,] %>% 
  dplyr::filter(is.na(pySCENIC_Regulon)) %>% select(TF) %>% deframe()
DotPlot(seuratobj[,seuratobj$ident == "Thymus"], features=x, group.by="Annotation") + 
  RotatedAxis() + scale_color_viridis() + 
  xlab(paste0("Integrative regulon specific Thymus TFs\n(Not found by pySCENIC)")) +
  theme(panel.background = element_rect(fill = "black"))
ggsave("/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_TF_expression_Thymus_specificIntegrative.pdf", width = 7, height = 3.5)



x <- regulonID_df[regulonID_df$TF %in% SSFintRegul$Liver,] %>% 
  dplyr::filter(is.na(pySCENIC_Regulon)) %>% select(TF) %>% deframe()
DotPlot(seuratobj[,seuratobj$ident == "Liver"], features=x, group.by="Annotation") + 
  RotatedAxis() + scale_color_viridis() + 
  xlab(paste0("Integrative regulon specific Liver TFs\n(Not found by pySCENIC)")) +
  theme(panel.background = element_rect(fill = "black"))
ggsave("/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_TF_expression_Liver_specificIntegrative.pdf", width = 7, height = 5)



#------------------------------------------------------------------------------#
#                     Expression of each gene by tissue                        #
#------------------------------------------------------------------------------#

allSSF_gg <-  lapply(names(SSFintRegul), function(id){
  ssf <- SSFintRegul[[id]]
  # Get TF identified but our method but not pySCENIC
  x <- regulonID_df[regulonID_df$TF %in% ssf,] %>% 
    select(TF) %>% deframe()
  
  p <- DotPlot(seuratobj, features=x, group.by="ident", cluster.idents = TRUE) + RotatedAxis() + 
    scale_color_viridis() + 
    xlab(paste0("Integrative regulon specific ", id, " TFs")) +
    theme(panel.background = element_rect(fill = "black"))
  return(p)
})
patchwork::wrap_plots(allSSF_gg)
ggsave("/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_TF_expression_tissue_specific.pdf", width = 17, height = 7)
