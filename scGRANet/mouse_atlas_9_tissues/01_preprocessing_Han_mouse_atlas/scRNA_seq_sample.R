library(tidyverse)

#------------------------------------------------------------------------------#
#                                 Create sparse matrix                         #
#------------------------------------------------------------------------------#
# It was not possible to read hd5 file into R
# the file was readed into python and:
# import scanpy
# x = scanpy.read_h5ad('/home/bq_aquintero/projects/mouse_brain_TF_activity/data/scrna/Han_atlas/MCA_BatchRemoved_Merge_dge.h5ad')
# x.write_csvs('/home/bq_aquintero/projects/mouse_brain_TF_activity/data/scrna/Han_atlas/MCA_BatchRemoved_Merge_dge.csv', skip_data = False)

# # Read data
# # Cell names
# x1 <- read_csv("~/projects/mouse_brain_TF_activity/data/scrna/Han_atlas/MCA_BatchRemoved_Merge_dge/obs.csv") # cell Names
# # Gene Names
# x2 <- read_csv("~/projects/mouse_brain_TF_activity/data/scrna/Han_atlas/MCA_BatchRemoved_Merge_dge/var.csv") # Gene names
# # data
# x <- read_csv("~/projects/mouse_brain_TF_activity/data/scrna/Han_atlas/MCA_BatchRemoved_Merge_dge/X.csv", col_names=x2$index) # Gene names
# x <- as.matrix(x)
# rownames(x) <- x1$index
#
# x_sparse <- Matrix(x, sparse = TRUE)
# saveRDS(x_sparse, "~/projects/mouse_brain_TF_activity/data/scrna/Han_atlas/MCA_BatchRemoved_Merge_dge.RDS")

#------------------------------------------------------------------------------#
#                              Save data 9 tissues                             #
#------------------------------------------------------------------------------#
x_sparse <- readRDS("~/projects/mouse_atlas_TF_activity/data/scrna/Han_atlas/MCA_BatchRemoved_Merge_dge.RDS")
x_sparse

# Read cell annotation
cellAnnot <- read_csv("~/projects/mouse_atlas_TF_activity/data/scrna/Han_atlas/MCA_BatchRemoved_Merge_dge_cellinfo.csv")

# Read cell annotation
cellID <- read_csv("~/projects/mouse_atlas_TF_activity/data/scrna/Han_atlas/MCA_CellAssignments.csv")

cellAnnot %>%
  rename(Cell.name = index) %>%
  left_join(cellID,  by="Cell.name")

cellAnnot <- cellAnnot %>%
  rename(Cell.name = index) %>%
  left_join(cellID,  by="Cell.name")


byTissue <- cellAnnot %>%
  group_by(tissue) %>%
  summarise(n = n())

tissuesList <- c(
  "AdultBrain",
  "AdultLung",
  "AdultKidney",
  "AdultBoneMarrow",
  "AdultLiver",
  "AdultThymus",
  "AdultSpleen",
  #"AdultTestis",
  "AdultSmallIntestine"
)
# Select only desired tissues
unique(cellAnnot$tissue)
table(is.na(cellAnnot$Annotation))

cellAnnot <- cellAnnot %>%
  dplyr::filter(tissue %in% tissuesList)
unique(cellAnnot$tissue)
table(is.na(cellAnnot$Annotation))

unique(cellAnnot$Annotation)

table(rownames(x_sparse) %in% cellAnnot$Cell.name)
table(cellAnnot$Cell.name %in% rownames(x_sparse))

idx <- match(cellAnnot$Cell.name, rownames(x_sparse))
x_sparse <- x_sparse[idx, ]
dim(x_sparse)

cellAnnot <- as.data.frame(cellAnnot)
rownames(cellAnnot) <- cellAnnot$Cell.name

#------------------------------------------------------------------------------#
#                             Create Seurat object.                            #
#------------------------------------------------------------------------------#
library(Seurat)
# Create seurat object
seuratobj <- CreateSeuratObject(counts = t(x_sparse),
                                min.cells = 3, min.features = 200,
                                meta.data = cellAnnot)
seuratobj
dim(seuratobj)
dim(seuratobj@meta.data)
head(seuratobj@meta.data)



#saveRDS(seuratobj, "~/projects/mouse_atlas_TF_activity/data/scrna/Han_atlas/MCA_Adult_9Tissues_Seurat.RDS")


seuratobj@meta.data %>%
  ggplot(aes(x = tissue)) +
  geom_bar()


#------------------------------------------------------------------------------#
#                   Filter genes expressed in few cells                        #
#------------------------------------------------------------------------------#
# isexpressed <- rowSums(seuratobj@assays$RNA@counts > 0)
# hist(isexpressed)
# table(isexpressed > 0)
# # expressed in at least 5 percent of the cells
# table(isexpressed > ncol(seuratobj)*0.10)
#

percent <- 0.01
minCountsPerGene <- 3 * percent * ncol(seuratobj@assays$RNA@counts)
nCountsPerGene <- rowSums(seuratobj@assays$RNA@counts, na.rm = T)
table(nCountsPerGene > minCountsPerGene)

minSamples <- ncol(seuratobj@assays$RNA@counts) * percent
nCellsPerGene <- rowSums(seuratobj@assays$RNA@counts>0, na.rm = T)
table(nCellsPerGene > minSamples)

table((nCellsPerGene > minSamples) & (nCountsPerGene > minCountsPerGene))
table((nCellsPerGene > minSamples) + (nCountsPerGene > minCountsPerGene))



idx <- (nCellsPerGene > minSamples) & (nCountsPerGene > minCountsPerGene)
table(idx)
seuratobj <- seuratobj[idx,]
seuratobj


# Check if all cells have expression
isexpressed <- colSums(seuratobj@assays$RNA@counts > 0)
hist(isexpressed)
table(isexpressed > 0)


#------------------------------------------------------------------------------#
#                               sample dataset                                 #
#------------------------------------------------------------------------------#
seuratobj_1K  <- seuratobj[,sort(sample(ncol(seuratobj), size=1000))]
seuratobj_5K  <- seuratobj[,sort(sample(ncol(seuratobj), size=5000))]
seuratobj_10K <- seuratobj[,sort(sample(ncol(seuratobj), size=10000))]
seuratobj_1K
seuratobj_5K
seuratobj_10K


out_path <- "/home/bq_aquintero/projects/GRANet/mouse_atlas_8tissues/data/scrna"
dir.create(out_path, recursive = TRUE)


saveRDS(seuratobj,     file.path(out_path, "mouse_scrna_42K.RDS"))
saveRDS(seuratobj_1K,  file.path(out_path, "mouse_scrna_1K.RDS"))
saveRDS(seuratobj_5K,  file.path(out_path, "mouse_scrna_5K.RDS"))
saveRDS(seuratobj_10K, file.path(out_path, "mouse_scrna_10K.RDS"))

seuratobj@meta.data %>%
  ggplot(aes(x = tissue)) +
  geom_bar() +
seuratobj_1K@meta.data %>%
  ggplot(aes(x = tissue)) +
  geom_bar() +
seuratobj_5K@meta.data %>%
  ggplot(aes(x = tissue)) +
  geom_bar() +
seuratobj_10K@meta.data %>%
  ggplot(aes(x = tissue)) +
  geom_bar()
