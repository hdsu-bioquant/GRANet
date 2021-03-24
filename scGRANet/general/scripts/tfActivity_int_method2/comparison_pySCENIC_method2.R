library(msigdbr)
library(clusterProfiler)
library(ggupset)

p_exprs              <- 'results/scrna/seurat/rna_seurat_int_annot_transfer.RDS'
p_nmfAUCpySCENIC <- "results/scrna/SCENIC/NMF_auc_mtx.RDS"
p_nmfregulonAUC  <- 'results/integrated/TF_activity_method2_5000/NMF_regulonAUC.RDS'
p_regulonpySCENIC <- "results/scrna/SCENIC/regulons.csv"
p_regulon <- "results/integrated/TF_activity_method2_5000/tfRegulons_asDF.RDS"

#------------------------------------------------------------------------------#
#                               Read Data                                      #
#------------------------------------------------------------------------------#
  
exprs <- readRDS(p_exprs)
nmfAUCpySCENIC <- readRDS(p_nmfAUCpySCENIC)
nmfregulonAUC  <- readRDS(p_nmfregulonAUC)


#------------------------------------------------------------------------------#
#                                Read regulons                                 #
#------------------------------------------------------------------------------#

regulonByCell <- readRDS(p_regulon)
regulonpySCENIC <- read_csv(p_regulonpySCENIC, col_names = TRUE, skip = 1) %>% 
  dplyr::rename(TF = X1) %>% 
  dplyr::slice(-1) %>% 
  dplyr::select(TF, TargetGenes) %>% 
  #head() %>% 
  mutate(TargetGenes = gsub("\\[|\\('|\\)|'|\\]| ", "", TargetGenes)) %>% 
  mutate(TargetGenes = strsplit(TargetGenes, ",")) %>% 
  mutate(TargetGenes = lapply(TargetGenes, function(x) x[seq(1, length(x), by=2)])) %>% 
  group_by(TF) %>% 
  summarise(TargetGenes = do.call(c, TargetGenes)) %>% 
  distinct()


regulonByCell

#------------------------------------------------------------------------------#
#                  Signature specific features & H heatmap                     #
#------------------------------------------------------------------------------#

ssf_pySCENIC <- SignatureSpecificFeatures(nmfAUCpySCENIC, k = 8,return_all_features = TRUE)
ssf_integrative <- SignatureSpecificFeatures(nmfregulonAUC, k = 8,return_all_features = TRUE)
rownames(ssf_pySCENIC) <- sub("\\(.*", "", rownames(ssf_pySCENIC))
rownames(ssf_integrative) <- sub(" \\(.*", "", rownames(ssf_integrative))


full_join(tibble(TF=rownames(ssf_pySCENIC),
                 pySCENIC=ssf_pySCENIC[,5]),
          
          tibble(TF=rownames(ssf_integrative ),
                 Integrative=ssf_integrative [,5])) %>% 
  #mutate(pySCENIC = pySCENIC+1) %>% 
  #mutate(Integrative = Integrative+1) %>% 
  mutate(pySCENIC = if_else(is.na(pySCENIC), 0, pySCENIC)) %>%
  mutate(Integrative = if_else(is.na(Integrative), 0, Integrative)) %>% 
  pivot_longer(cols = -TF, names_to = "SigID", values_to = "IsSig") %>% 
  dplyr::filter(IsSig == 1 ) %>% 
  dplyr::select(-IsSig ) %>% 
  group_by(TF) %>%
  summarize(SigID = list(SigID)) %>% 
  ggplot(aes(x = SigID)) +
  #geom_bar(fill=c(rep("red",8),rep("black",12))) +
  geom_bar() +
  ggtitle("Neutrophil defining regulons") +
  scale_x_upset(order_by = "degree", n_intersections = 3) +
  cowplot::theme_cowplot()
  

# Ciliated cells
# py 3
# int 7
full_join(tibble(TF=rownames(ssf_pySCENIC),
                 pySCENIC=ssf_pySCENIC[,3]),
          
          tibble(TF=rownames(ssf_integrative ),
                 Integrative=ssf_integrative [,7])) %>% 
  #mutate(pySCENIC = pySCENIC+1) %>% 
  #mutate(Integrative = Integrative+1) %>% 
  mutate(pySCENIC = if_else(is.na(pySCENIC), 0, pySCENIC)) %>%
  mutate(Integrative = if_else(is.na(Integrative), 0, Integrative)) %>% 
  pivot_longer(cols = -TF, names_to = "SigID", values_to = "IsSig") %>% 
  dplyr::filter(IsSig == 1 ) %>% 
  dplyr::select(-IsSig ) %>% 
  group_by(TF) %>%
  summarize(SigID = list(SigID)) %>% 
  ggplot(aes(x = SigID)) +
  #geom_bar(fill=c(rep("red",8),rep("black",12))) +
  geom_bar() +
  ggtitle("Ciliated defining regulons") +
  scale_x_upset(order_by = "degree", n_intersections = 3) +
  cowplot::theme_cowplot() 


ssf_gg <- ss_features %>% 
  as_tibble(rownames = "geneID") %>% 
  pivot_longer(cols = -geneID, names_to = "SigID", values_to = "IsSig") %>% 
  dplyr::filter(IsSig == 1 ) %>% 
  dplyr::select(-IsSig ) %>% 
  group_by(geneID) %>%
  summarize(SigID = list(SigID)) %>% 
  ggplot(aes(x = SigID)) +
  geom_bar(fill=c(rep("red",8),rep("black",12))) +
  scale_x_upset(order_by = "degree", n_intersections = 20) +
  cowplot::theme_cowplot()
ssf_gg


# CTL cells
# py 6
# int 2
full_join(tibble(TF=rownames(ssf_pySCENIC),
                 pySCENIC=ssf_pySCENIC[,6]),
          
          tibble(TF=rownames(ssf_integrative ),
                 Integrative=ssf_integrative [,2])) %>% 
  #mutate(pySCENIC = pySCENIC+1) %>% 
  #mutate(Integrative = Integrative+1) %>% 
  mutate(pySCENIC = if_else(is.na(pySCENIC), 0, pySCENIC)) %>%
  mutate(Integrative = if_else(is.na(Integrative), 0, Integrative)) %>% 
  pivot_longer(cols = -TF, names_to = "SigID", values_to = "IsSig") %>% 
  dplyr::filter(IsSig == 1 ) %>% 
  dplyr::select(-IsSig ) %>% 
  group_by(TF) %>%
  summarize(SigID = list(SigID)) %>% 
  ggplot(aes(x = SigID)) +
  #geom_bar(fill=c(rep("red",8),rep("black",12))) +
  geom_bar() +
  ggtitle("CTL defining regulons") +
  scale_x_upset(order_by = "degree", n_intersections = 3) +
  cowplot::theme_cowplot() 


ssf_gg <- ss_features %>% 
  as_tibble(rownames = "geneID") %>% 
  pivot_longer(cols = -geneID, names_to = "SigID", values_to = "IsSig") %>% 
  dplyr::filter(IsSig == 1 ) %>% 
  dplyr::select(-IsSig ) %>% 
  group_by(geneID) %>%
  summarize(SigID = list(SigID)) %>% 
  ggplot(aes(x = SigID)) +
  geom_bar(fill=c(rep("red",8),rep("black",12))) +
  scale_x_upset(order_by = "degree", n_intersections = 20) +
  cowplot::theme_cowplot()
ssf_gg



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
          show_row_names = TRUE,
          show_column_dend = FALSE,
          show_row_dend = FALSE,
          cluster_rows = TRUE)
}


hAUCpySCENIC <- HMatrix(nmfAUCpySCENIC, k=8)
Heatmap_Regulon(hAUCpySCENIC, exprs, TRUE, "H AUCpySCENIC")
hregulonAUC <- HMatrix(nmfregulonAUC, k=8)
rownames(hregulonAUC) <- paste0("Sign.", 1:nrow(hregulonAUC))
Heatmap_Regulon(hregulonAUC, exprs, TRUE, "H regulonAUC")


#------------------------------------------------------------------------------#
#                        signature expuse UMAP.                                #
#------------------------------------------------------------------------------#


Seurat::FetchData(exprs, vars=c("predicted.id", "UMAP_1", "UMAP_2")) %>% 
  rownames_to_column("CellID") %>% 
  right_join(rownames_to_column(as.data.frame(t(hregulonAUC)), "CellID")) %>% 
  pivot_longer(cols = -c("UMAP_1", "UMAP_2", "CellID", "predicted.id"), names_to = "SigID", values_to = "Exposure") %>% 
  arrange(Exposure) %>% 
  # Normalize by signature
  group_by(SigID) %>% 
  mutate(Exposure = Exposure/max(Exposure)) %>% 
  
  ggplot(aes(x=UMAP_1, y=UMAP_2, color=Exposure)) +
  geom_point() +
  scale_color_viridis() +
  facet_wrap(.~SigID) +
  cowplot::theme_cowplot()




Seurat::FetchData(exprs, vars=c("predicted.id", "UMAP_1", "UMAP_2")) %>% 
  rownames_to_column("CellID") %>% 
  right_join(rownames_to_column(as.data.frame(t(hregulonAUC)), "CellID")) %>% 
  pivot_longer(cols = -c("UMAP_1", "UMAP_2", "CellID", "predicted.id"), names_to = "SigID", values_to = "Exposure") %>% 
  arrange(Exposure) %>% 
  # Normalize by signature
  group_by(SigID) %>% 
  mutate(Exposure = Exposure/max(Exposure)) %>% 
  # Keep only Secretory Signatures
  dplyr::filter(SigID %in% c("Sign.3", "Sign.4", "Sign.6")) %>% 
  
  ggplot(aes(x=UMAP_1, y=UMAP_2, color=Exposure)) +
  geom_point() +
  scale_color_viridis() +
  facet_wrap(.~SigID) +
  cowplot::theme_cowplot()




Seurat::GetAssayData(exprs, slot="scale.data")
Seurat::GetAssay(exprs)

exprs@meta.data



#------------------------------------------------------------------------------#
#                       Top Q features by signature                            #
#------------------------------------------------------------------------------#


#Neutrophilpyscenic 5
#Neutrophil method2 5

top_Qperc_assing <- function(wmatrix, Q=0.8){
  colnames(wmatrix) <- paste0("Sign.", 1:ncol(wmatrix))
  sig_assign <- lapply(setNames(colnames(wmatrix), colnames(wmatrix)), function(sigID){
    selec_wmatrix <- do.call(cbind, lapply(as.data.frame(wmatrix), function(sign_expo){
      sign_expo[sign_expo < quantile(sign_expo, Q)] <- NA
      sign_expo
    }))
    rownames(selec_wmatrix) <- rownames(wmatrix)
    selec_wmatrix <- selec_wmatrix[!is.na(selec_wmatrix[,sigID]),,drop=FALSE]
    # Keep only the top feature if there's an overlap
    sig_SE_IDs <- rownames(selec_wmatrix[rowMaxs(selec_wmatrix, na.rm = TRUE) == selec_wmatrix[,sigID],])
    sig_SE_IDs
  })
  sig_assign
}


wAUCpySCENIC <- WMatrix(nmfAUCpySCENIC, k=8)
wregulonAUC <- WMatrix(nmfregulonAUC, k=8)
rownames(wAUCpySCENIC) <- sub("\\(.*", "", rownames(wAUCpySCENIC))
rownames(wregulonAUC) <- sub(" \\(.*", "", rownames(wregulonAUC))

sf_AUCpySCENIC <- top_Qperc_assing(wAUCpySCENIC, 0.7)
sf_regulonAUC <- top_Qperc_assing(wregulonAUC, 0.7)




#------------------------------------------------------------------------------#
#                    Neutrophile signature term erichment                      #
#------------------------------------------------------------------------------#


msigdb_hs <- msigdbr(species = "Homo sapiens")
term2gene <- msigdb_hs %>% 
  mutate(term = gs_name) %>% 
  mutate(gene = gene_symbol) %>% 
  dplyr::select(term, gene)


neutrophil_regulonTop <- list(pySCENIC = unique(regulonpySCENIC$TargetGenes[regulonpySCENIC$TF %in% sf_AUCpySCENIC$Sign.5]),
                              Integrative = unique(regulonByCell$Neu$target[regulonByCell$Neu$TF %in% sf_regulonAUC$Sign.5]))

neutrophil_regulonTop <- list(pySCENIC = unique(sf_AUCpySCENIC$Sign.5),
                              Integrative = unique(sf_regulonAUC$Sign.5))



idx <- grep(pattern="imm", term2gene$term, ignore.case=TRUE) 
term2gene <- term2gene[idx,]

sign_compare_t10_Msig <- compareCluster(geneClusters = neutrophil_regulonTop,
                                        fun = "enricher",
                                        TERM2GENE = term2gene)
dotplot(sign_compare_t10_Msig, showCategory = 10)


# sign_compare_t10_Msig <- compareCluster(geneClusters = sign_features, 
#                                         fun = "enricher",
#                                         TERM2GENE = term2gene)
# dotplot(sign_compare_t10_Msig, showCategory = 5)


list(pySCENIC = sf_AUCpySCENIC$`Sign. 5`,
Integrative = sf_regulonAUC$`Sign. 5`)




#------------------------------------------------------------------------------#
#                    Compare W Neutrophile signature                           #
#------------------------------------------------------------------------------#

full_join(tibble(TF=rownames(wAUCpySCENIC),
                 pySCENIC=wAUCpySCENIC[,5]),
          
          tibble(TF=rownames(wregulonAUC ),
                 Integrative=wregulonAUC [,5])) %>% 
  #mutate(pySCENIC = pySCENIC+1) %>% 
  #mutate(Integrative = Integrative+1) %>% 
  mutate(pySCENIC = if_else(is.na(pySCENIC), 0, pySCENIC)) %>%
  mutate(Integrative = if_else(is.na(Integrative), 0, Integrative)) %>% 

  ggplot(aes(x=pySCENIC, y=Integrative)) +
    geom_point()
          




