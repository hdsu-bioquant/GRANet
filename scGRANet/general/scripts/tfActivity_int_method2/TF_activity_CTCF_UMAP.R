options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

p_exprs              <- as.character(args[1])
p_regulonAUC_devCor  <- as.character(args[2])
p_CTCF_TFactivity    <- as.character(args[3])

#------------------------------------------------------------------------------#
p_exprs              <- 'results/scrna/seurat/rna_seurat_int_annot_transfer.RDS'
p_regulonAUC_devCor  <- "results/integrated/TF_activity_method2_500000/regulonAUC_devCorrected.RDS"
p_CTCF_TFactivity    <- "figures/tfActivity_integrated/method2_500000_TFactivity_CTCF_devCorrected.pdf"
#------------------------------------------------------------------------------#

library(tidyverse)
library(Seurat)
library(patchwork)

##-----------------------------------------------------------------------------#
##                        Regulon heatmap annotation                           #
##-----------------------------------------------------------------------------#
exprs         <- readRDS(p_exprs)
TFactivity    <- readRDS(p_regulonAUC_devCor)
#------------------------------------------------------------------------------#
#                             Heatmap pyScenic regulon AUC                     #
#------------------------------------------------------------------------------#

grep("CTCF", rownames(TFactivity), value = TRUE)

exprs_fil <- exprs[,colnames(exprs) %in% colnames(TFactivity)]


p1 <- as.data.frame(exprs_fil@reductions$umap@cell.embeddings) %>%
  rownames_to_column("CellID") %>%
  mutate(SelTF = TFactivity["CTCF (+)", match(CellID, colnames(TFactivity))]) %>%
  mutate(SelTF = SelTF/max(SelTF)) %>%
  arrange(SelTF) %>% 
  ggplot(aes(x=UMAP_1, y=UMAP_2, color=SelTF)) +
  geom_point() +
  scale_color_viridis() +
  labs(color="CTCF (+)") +
  cowplot::theme_cowplot()

p2 <- as.data.frame(exprs_fil@reductions$umap@cell.embeddings) %>%
  rownames_to_column("CellID") %>%
  mutate(CellType = exprs_fil$predicted.id) %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, color=CellType)) +
  geom_point() +
  cowplot::theme_cowplot()


p <- p1+p2
ggsave(filename = p_CTCF_TFactivity, plot = p, 
       width = 9, height = 3.5)

