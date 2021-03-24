
geneAnnot <- read_csv("~/projects/mouse_brain_TF_activity/data/scrna/gene_annotate.csv")
table(duplicated(geneAnnot$gene_short_name))

cellAnnot <- read_csv("~/projects/mouse_brain_TF_activity/data/scrna/cell_annotate.csv")
head(cellAnnot)
unique(cellAnnot$Main_cell_type)
#unique(cellAnnot$)

scrna <- readRDS("~/projects/mouse_brain_TF_activity/data/scrna/gene_count_cleaned_sampled_100k.RDS")
scrna


dim(scrna)
