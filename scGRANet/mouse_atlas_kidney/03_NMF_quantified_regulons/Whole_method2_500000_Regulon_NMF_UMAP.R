params <- list(
  K                   = 11,
  p_exprs             = '/home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scrna/Han_atlas/MCA_Adult_9Tissues_Seurat.RDS',
  p_nmf = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_regulonAUC.RDS",
  p_umap = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_regulonAUC_UMAP.pdf",
  annotCol = "tissue",
  p_archrproj = "~/projects/mouse_atlas_TF_activity/adult_9tissues/results/scatac/archr/ArchR02_MotifMatch/")


norm_nmf_exp <- readRDS(params$p_nmf)
exprs <- readRDS(params$p_exprs)

colID <- params$annotCol
type_colVector <- exprs@meta.data[,colID, drop=FALSE] %>% 
  dplyr::distinct() %>% 
  #dplyr::rename(Celltype =  !!enquo(colID)) %>% 
  dplyr::rename(Celltype =  !!colID) %>% 
  dplyr::arrange(Celltype) %>% 
  dplyr::mutate(color = alphabet.colors(n())) %>% 
  deframe()

##----------------------------------------------------------------------------##
##                         UMAP H matrix                                      ##
##----------------------------------------------------------------------------##
hmatrix_norm <- HMatrix(norm_nmf_exp, k = params$K)
umapView <- umap(t(hmatrix_norm))

umapView_df <- as.data.frame(umapView$layout)
colnames(umapView_df) <- c("UMAP1", "UMAP2")

p1 <- umapView_df %>% 
  rownames_to_column("cellID") %>% 
  left_join(rownames_to_column(exprs@meta.data, "cellID"), by = "cellID") %>% 
  ggplot(aes(x=UMAP1, y=UMAP2, color = !!sym(params$annotCol))) + 
  geom_point(size = .1, alpha = 0.95) + 
  scale_color_manual(values = type_colVector) +
  cowplot::theme_cowplot()
p1




##----------------------------------------------------------------------------##
##                         Add archr UMAP                                     ##
##----------------------------------------------------------------------------##


# Ncores      <- 20
#------------------------------------------------------------------------------#
library(ArchR)
library(patchwork)
library(ComplexHeatmap)
library(viridis)

ArchR::addArchRThreads(Ncores)
#addArchRGenome(genome)


#------------------------------------------------------------------------------#
#                          Debugged functions                                  #
#------------------------------------------------------------------------------#

archrproj <- loadArchRProject(path = params$p_archrproj)

atacdf <- archrproj@embeddings$UMAP$df
colnames(atacdf) <- c("UMAP1", "UMAP2")
atacdf$tissue <- archrproj$predictedGroup_Un

p2 <- atacdf %>% 
  ggplot(aes(x=UMAP1, y=UMAP2, color = tissue)) + 
  geom_point(size = .1, alpha = 0.95) + 
  scale_color_manual(values = type_colVector) +
  cowplot::theme_cowplot()
p2

p <- p1 + p2

ggsave(params$p_umap, p, width = 20, height = 10)
