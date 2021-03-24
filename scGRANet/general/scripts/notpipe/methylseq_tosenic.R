seuratobj <- readRDS('/media/ag-cherrmann/cherrmann/MethylSeq.rds')
dim(seuratobj)
exprs2scenic <- seuratobj@assays$RNA@counts
dim(exprs2scenic)
exprs2scenic <- t(exprs2scenic)
dim(exprs2scenic)
write.table(exprs2scenic, file = '/media/ag-cherrmann/cherrmann/MethylSeq_counts.tsv', 
            quote = FALSE, sep = '\t', row.names = TRUE, col.names = TRUE)





covid_nbt_main <- NormalizeData(covid_nbt_main, normalization.method = "LogNormalize", scale.factor = 10000)
covid_nbt_main <- ScaleData(covid_nbt_main, features = rownames(covid_nbt_main))
covid_nbt_main <- FindVariableFeatures(covid_nbt_main, selection.method = "vst", nfeatures = 2000)





seuratobj <- NormalizeData(seuratobj, normalization.method = "LogNormalize", scale.factor = 10000)
seuratobj <- ScaleData(seuratobj, features = rownames(seuratobj))
dim(seuratobj)
seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst", nfeatures = 8000)

var_genes <- VariableFeatures(seuratobj)
exprs2scenic <- GetAssayData(seuratobj)[var_genes,]
dim(exprs2scenic)

exprs2scenic <- t(exprs2scenic)
dim(exprs2scenic)
write.table(exprs2scenic, file = '/media/ag-cherrmann/cherrmann/MethylSeq_SCENIC/MethylSeq_counts_8000hvg.tsv', 
            quote = FALSE, sep = '\t', row.names = TRUE, col.names = TRUE)


##-----------------------------------------------------------------------------#
##                        Regulon 8000 HVG                                     #
##-----------------------------------------------------------------------------#


# SCENIC commands:
#singularity shell --bind /home/bq_aquintero/10_charite_covid19/ --bind /media/ag-cherrmann/cherrmann/MethylSeq_SCENIC/  docker://aertslab/pyscenic:0.10.4
#cd /media/ag-cherrmann/cherrmann/MethylSeq_SCENIC/
#arboreto_with_multiprocessing.py -o expr_mat.adjacencies.tsv --num_workers 20 MethylSeq_counts_8000hvg.tsv /home/bq_aquintero/10_charite_covid19/src/charite_covid19/aux/scenic/allTFs_hg38.txt

# pyscenic ctx expr_mat.adjacencies.tsv /home/bq_aquintero/10_charite_covid19/src/charite_covid19/aux/scenic/rcistarget/hg19-tss-centered-10kb-7species.mc9nr.feather  /home/bq_aquintero/10_charite_covid19/src/charite_covid19/aux/scenic/rcistarget/hg19-tss-centered-5kb-7species.mc9nr.feather \
# --annotations_fname /home/bq_aquintero/10_charite_covid19/src/charite_covid19/aux/scenic/rcistarget/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
# --expression_mtx_fname MethylSeq_counts_8000hvg.tsv \
# --mode "dask_multiprocessing" \
# --output regulons.csv \
# --num_workers 20

# pyscenic aucell MethylSeq_counts_8000hvg.tsv regulons.csv \
# -o auc_mtx.csv \
# --num_workers 20

##-----------------------------------------------------------------------------#
##                        Regulon heatmap annotation                           #
##-----------------------------------------------------------------------------#


regulonAUC_pySCENIC <- read_csv('/media/ag-cherrmann/cherrmann/MethylSeq_SCENIC/scenic_8000hvg/auc_mtx.csv') %>% 
  column_to_rownames("Cell") %>% 
  t()
dim(regulonAUC_pySCENIC)

unique(seuratobj$celltype_proposal)


Heatmap_Regulon <- function(regulon, seurat, normRows, title){
  
  type.colVector <- seurat@meta.data[,"celltype_proposal", drop=FALSE] %>% 
    dplyr::distinct() %>% 
    dplyr::rename(Celltype =  celltype_proposal) %>% 
    dplyr::arrange(Celltype) %>% 
    dplyr::mutate(color = alphabet.colors(n())) %>% 
    deframe()
  
  # Build Heatmap annotation
  heat.anno <- HeatmapAnnotation(df  = seurat@meta.data[colnames(regulon),"celltype_proposal", drop=FALSE],
                                 col = list(celltype_proposal = type.colVector),
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
          show_row_names = FALSE,
          show_column_dend = FALSE,
          show_row_dend = FALSE,
          use_raster = TRUE,
          cluster_rows = TRUE)
}

p <- Heatmap_Regulon(regulon=regulonAUC_pySCENIC, seurat=seuratobj, normRows=TRUE,
                     title = "Regulon AUC\npySCENIC\n8000 HVG")

pdf('/media/ag-cherrmann/cherrmann/MethylSeq_SCENIC/scenic_8000hvg/heatmap_auc_mtx.pdf', width = 10, height = 6)
p
dev.off()


##-----------------------------------------------------------------------------#
##                        Regulon all genes.                                   #
##-----------------------------------------------------------------------------#


# SCENIC commands:

# Interactive session:
"qsub -I -l 'walltime=72:00:00, nodes=1:ppn=15,  mem=500g'"

#singularity shell --bind /home/bq_aquintero/10_charite_covid19/ --bind /media/ag-cherrmann/cherrmann/MethylSeq_SCENIC/  docker://aertslab/pyscenic:0.10.4
#cd /media/ag-cherrmann/cherrmann/MethylSeq_SCENIC/

# Step 1
#arboreto_with_multiprocessing.py -o expr_mat.adjacencies_allgenes.tsv --num_workers 15 MethylSeq_counts.tsv /home/bq_aquintero/10_charite_covid19/src/charite_covid19/aux/scenic/allTFs_hg38.txt

# Step 2
# pyscenic ctx expr_mat.adjacencies_allgenes.tsv /home/bq_aquintero/10_charite_covid19/src/charite_covid19/aux/scenic/rcistarget/hg19-tss-centered-10kb-7species.mc9nr.feather  /home/bq_aquintero/10_charite_covid19/src/charite_covid19/aux/scenic/rcistarget/hg19-tss-centered-5kb-7species.mc9nr.feather \
# --annotations_fname /home/bq_aquintero/10_charite_covid19/src/charite_covid19/aux/scenic/rcistarget/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
# --expression_mtx_fname MethylSeq_counts.tsv \
# --mode "dask_multiprocessing" \
# --output regulons.csv \
# --num_workers 15

# pyscenic aucell MethylSeq_counts.tsv regulons.csv \
# -o auc_mtx.csv \
# --num_workers 15

