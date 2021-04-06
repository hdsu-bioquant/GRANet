library(tidyverse)

#------------------------------------------------------------------------------#
#                            Read Kidney expression data                       #
#------------------------------------------------------------------------------#
data_path <- "/media/ag-cherrmann/projects/10_charite_covid19/subprojects/human_atlas_TF_activity/data/rna/"
data_loc <- tibble(path = list.files(file.path(data_path, "dge_rmbatch"), full.names = TRUE))
data_loc <- data_loc %>% 
  mutate(file = basename(path)) %>% 
  mutate(isKidney = grepl("lung", file, ignore.case = TRUE)) %>% 
  dplyr::filter(isKidney & grepl("adult", file, ignore.case = TRUE))

dge_list <- apply(data_loc, 1, function(x){
  #print(x)
  read.table(x["path"], header = TRUE)
})
lapply(dge_list, dim)
common_genes <- Reduce(intersect, lapply(dge_list, rownames))
dge_list <- lapply(dge_list, function(x){
  x[common_genes,]
})
lapply(dge_list, dim)
dge <- do.call("cbind", dge_list)
dim(dge)
dge[1:5,1:5]

#------------------------------------------------------------------------------#
#                            Read Kidney annotation data                       #
#------------------------------------------------------------------------------#
data_loc <- tibble(path = list.files(file.path(data_path, "annotation_rmbatch_data_revised417"), full.names = TRUE))
data_loc <- data_loc %>% 
  mutate(file = basename(path)) %>% 
  mutate(isTissue = grepl("lung", file, ignore.case = TRUE)) %>% 
  dplyr::filter(isTissue & grepl("adult", file, ignore.case = TRUE))

annot_list <- apply(data_loc, 1, function(x){
  read_csv(x["path"])
})
annot <- bind_rows(annot_list)
annot <- as.data.frame(annot)
rownames(annot) <- annot$Cell_id
table(colnames(dge) %in% annot$Cell_id)


#------------------------------------------------------------------------------#
#                            Create Seurat object.                             #
#------------------------------------------------------------------------------#
library(Seurat)
seuratobj <- CreateSeuratObject(dge[,annot$Cell_id], meta.data = annot)

rownames(seuratobj)
dir.create("~/projects/GRANet/human_atlas_lung/data/scrna/", recursive = TRUE)
saveRDS(seuratobj, "~/projects/GRANet/human_atlas_lung/data/scrna/adult_human_lung_seurat.RDS")


#------------------------------------------------------------------------------#
#                            Save data for pySCENIC                            #
#------------------------------------------------------------------------------#
exprs2scenic <- seuratobj@assays$RNA@counts
dim(exprs2scenic)

exprs2scenic <- t(as.matrix(exprs2scenic))
dim(exprs2scenic)
dir.create("~/projects/GRANet/human_atlas_lung/results/scrna/SCENIC/", recursive = TRUE)
write.table(exprs2scenic, file = "~/projects/GRANet/human_atlas_lung/results/scrna/SCENIC/rnaseq_counts.tsv", 
            quote = FALSE, sep = '\t', row.names = TRUE, col.names = TRUE)
