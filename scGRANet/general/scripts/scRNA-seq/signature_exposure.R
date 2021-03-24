library(msigdbr)
library(clusterProfiler)
library(ggupset)

p_exprs              <- 'results/scrna/seurat/rna_seurat_int_annot_transfer.RDS'
p_nmfexprs  <- 'results/scrna/seurat/NMF_geneExpresion_5000.RDS'
#p_regulonpySCENIC <- "results/scrna/SCENIC/regulons.csv"
#p_regulon <- "results/integrated/TF_activity_method2_5000/tfRegulons_asDF.RDS"

#------------------------------------------------------------------------------#
#                               Read Data                                      #
#------------------------------------------------------------------------------#

exprs <- readRDS(p_exprs)
nmfexprs <- readRDS(p_nmfexprs)
#nmfregulonAUC  <- readRDS(p_nmfregulonAUC)


#------------------------------------------------------------------------------#
#                               H heatmap                                      #
#------------------------------------------------------------------------------#
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


hnmfexprs <- HMatrix(nmfexprs, k=8)
rownames(hnmfexprs) <- paste0("Sign.", 1:nrow(hnmfexprs))
Heatmap_Regulon(hnmfexprs, exprs, TRUE, "H expression")


#------------------------------------------------------------------------------#
#                        signature expuse UMAP.                                #
#------------------------------------------------------------------------------#


Seurat::FetchData(exprs, vars=c("predicted.id", "UMAP_1", "UMAP_2")) %>% 
  rownames_to_column("CellID") %>% 
  right_join(rownames_to_column(as.data.frame(t(hnmfexprs)), "CellID")) %>% 
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
  right_join(rownames_to_column(as.data.frame(t(hnmfexprs)), "CellID")) %>% 
  pivot_longer(cols = -c("UMAP_1", "UMAP_2", "CellID", "predicted.id"), names_to = "SigID", values_to = "Exposure") %>% 
  arrange(Exposure) %>% 
  # Normalize by signature
  group_by(SigID) %>% 
  mutate(Exposure = Exposure/max(Exposure)) %>% 
  # Keep only Secretory Signatures
  dplyr::filter(SigID %in% c("Sign.1", "Sign.4", "Sign.6")) %>% 
  
  ggplot(aes(x=UMAP_1, y=UMAP_2, color=Exposure)) +
  geom_point() +
  scale_color_viridis() +
  facet_wrap(.~SigID) +
  cowplot::theme_cowplot()




