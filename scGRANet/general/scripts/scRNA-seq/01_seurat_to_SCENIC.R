options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

p_seuratobj    <- as.character(args[1])
p_countsSCENIC <- as.character(args[2])
hvg            <- tolower(as.character(args[3]))
#------------------------------------------------------------------------------#
# p_seuratobj    <-  "~/projects/mouse_atlas_TF_activity/data/scrna/Han_atlas/MCA_AdultLung_Seurat.RDS"
# p_countsSCENIC <-  "~/projects/mouse_atlas_TF_activity/adult_lung/results/scrna/SCENIC/rnaseq_counts.tsv"
# hvg  <- "all"
#------------------------------------------------------------------------------#
library(Seurat)

#------------------------------------------------------------------------------#
#                          Debugged functions                                  #
#------------------------------------------------------------------------------#
seuratobj <- readRDS(p_seuratobj)
dim(seuratobj@assays$RNA@counts)
#dim(seuratobj@assays$integrated@data)

if (!hvg == "all") {
  hvg <- as.numeric(hvg)
  if (nrow(seuratobj) <= hvg) {
    hvg  <- "all"
  }
}



if (hvg == "all") {
  #exprs2scenic <- GetAssayData(seuratobj)
  exprs2scenic <- seuratobj@assays$RNA@counts
  
} else {
  seuratobj <- NormalizeData(seuratobj, normalization.method = "LogNormalize", scale.factor = 10000)
  seuratobj <- ScaleData(seuratobj, features = rownames(seuratobj))
  seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst", nfeatures = hvg)
  var_genes <- VariableFeatures(seuratobj)
  exprs2scenic <- GetAssayData(seuratobj)[var_genes,]
  
}

dim(exprs2scenic)
exprs2scenic <- t(as.matrix(exprs2scenic))
dim(exprs2scenic)
write.table(exprs2scenic, file = p_countsSCENIC, 
            quote = FALSE, sep = '\t', row.names = TRUE, col.names = TRUE)
