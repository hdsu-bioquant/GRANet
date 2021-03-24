DATAPATH <- "/home/bq_aquintero/projects/charite_covid19_TF_activity/atac17_rna21/"
params <- list(regulonAUC = file.path(DATAPATH, 'results/integrated/TF_activity_method2_500000/regulonAUC.RDS'),
               seuratobj = "/home/bq_aquintero/projects/charite_covid19_TF_activity/results/scrna/seurat/rna_seurat_int_annot_transfer.RDS",
               p_archrproj = "/home/bq_aquintero/projects/charite_covid19_TF_activity/atac17_rna21/results/scatac/archr/ArchR02_MotifMatch/")

library(Seurat)
library(Polychrome)

reticulate::use_condaenv("tensor240pip", required = TRUE)
##----------------------------------------------------------------------------##
##                             read data.                                     ##
##----------------------------------------------------------------------------##


samples_df <- readxl::read_xlsx("/media/ag-cherrmann/projects/10_charite_covid19/data/scrna_seq_annot_path.xlsx") %>% 
  mutate(months_after_acute_phase = paste0(sub(" months", "", months_after_acute_phase), " months"))

regulon <- readRDS(params$regulonAUC)
regulon[1:5,1:5]
dim(regulon)

seuratobj <- readRDS(params$seuratobj)

head(seuratobj@meta.data)

idx <- match(seuratobj@meta.data$orig.ident, samples_df$SampleID)
seuratobj$Gender <- samples_df$Gender[idx]
seuratobj$Age <- samples_df$Age[idx]
seuratobj$months_after_acute_phase <- samples_df$months_after_acute_phase[idx]

##----------------------------------------------------------------------------##
##                             cell type by month                             ##
##----------------------------------------------------------------------------##

seuratobj@meta.data %>% 
  ggplot(aes(y= BioClassification)) +
  geom_bar()+
  facet_grid(.~months_after_acute_phase, scales="free") +
  ggtitle("Months after acute phase") +
  cowplot::theme_cowplot()


##----------------------------------------------------------------------------##
##                             regulon heatmap.                               ##
##----------------------------------------------------------------------------##

Heatmap_Regulon <- function(regulon, seurat, title){
  
  set.seed(1)
  idx <- sample(colnames(regulon), size=10000, replace=FALSE)
  regulon <- regulon[,idx]
  
  #colIDs <- c("Gender", "Age", "months_after_acute_phase")
  colIDs <- c("BioClassification", "Age", "months_after_acute_phase")
  names(colIDs) <- colIDs
  collist <- lapply(colIDs, function(colID){
    type.colVector <- seurat@meta.data[,colID, drop=FALSE] %>% 
      dplyr::distinct() %>% 
      #dplyr::rename(Celltype =  !!enquo(colID)) %>% 
      dplyr::rename(Celltype =  !!colID) %>% 
      dplyr::arrange(Celltype) %>% 
      dplyr::mutate(color = alphabet.colors(n())) %>% 
      deframe()
  })
  collist$months_after_acute_phase <- setNames(c("white", "yellow", "black"), c("0 months", "3 months", "6 months"))
  collist$Gender <- setNames(c("grey20", "grey80"), c("f", "m"))
  #circlize::colorRamp2(breaks = seq(from=min(seurat@meta.data$Age), to=max(seurat@meta.data$Age), length.out=10), colors = c("white", "black"))
  collist$Age <- circlize::colorRamp2(breaks = c(min(seurat@meta.data$Age), max(seurat@meta.data$Age)), colors = c("white", "black"))
  
  # 
  # return(collist)
  # seuratobj$Gender
  # seuratobj$Age
  # seuratobj$months_after_acute_phase
  
  
  
  # type.colVector <- seurat@meta.data[,colID, drop=FALSE] %>% 
  #   dplyr::distinct() %>% 
  #   #dplyr::rename(Celltype =  !!enquo(colID)) %>% 
  #   dplyr::rename(Celltype =  !!colID) %>% 
  #   dplyr::arrange(Celltype) %>% 
  #   dplyr::mutate(color = alphabet.colors(n())) %>% 
  #   deframe()
  # collist <- list(type.colVector) 
  # names(collist) <- colID
  
  # Build Heatmap annotation
  heat.anno <- HeatmapAnnotation(df  = seurat@meta.data[colnames(regulon),names(collist), drop=FALSE],
                                 col = collist,
                                 #list(predicted.id = type.colVector),
                                 show_annotation_name = TRUE, na_col = "white")
  
  
  
  
  
  #rownames(regulon) <- paste0("Sign.", 1:nrow(regulon))
  Heatmap(regulon,
          col = viridis(100),
          name = title,
          top_annotation = heat.anno,
          use_raster = TRUE,
          raster_device = "png",
          show_column_names = FALSE,
          show_row_names = FALSE,
          show_column_dend = FALSE,
          show_row_dend = FALSE,
          cluster_rows = TRUE)
}

h <- Heatmap_Regulon(regulon=regulon, seurat=seuratobj, title="Regulon AUC")

##----------------------------------------------------------------------------##
##                             NMF by cell type                               ##
##----------------------------------------------------------------------------##

celTypes_freq_df <- seuratobj@meta.data %>% 
  group_by(BioClassification) %>% 
  summarise(n = n()) %>% 
  dplyr::filter(n > 10)

lapply(celTypes_freq_df$BioClassification, function(id){
  
  seuratobjByID <- seuratobj[,seuratobj$BioClassification == id]
  print(table(seuratobjByID$BioClassification))
  
  regulonByID <- regulon[,colnames(regulon) %in% colnames(seuratobjByID)]
  #regulonByID[1:5,1:5]
  print(dim(regulonByID))
  
  if (ncol(regulonByID)>0) {
    
    outdir <- file.path(DATAPATH, "results/integrated/TF_activity_method2_500000", gsub(" ", "_", id))
    dir.create(outdir)
    print(outdir)
    #print(file.path(outdir, "regulonAUC.RDS"))
    
    rmarkdown::render(input = "src/charite_covid19/scripts/covid19/Regulon_NMF.Rmd",
                      output_file = file.path(outdir, paste0(gsub(" ", "_", id), "_regulonAUC.html")),
                       params = list(K     =  4,
                                     p_nmf = file.path(outdir, paste0("NMF_", gsub(" ", "_", id), "_regulonAUC.RDS"))
                                     ))
    
    
    # file.path(DATAPATH, "results/integrated/TF_activity_method2_500000/regulonAUC.RDS")
  }
  
  #sssss
  
  
})


dim(seuratobj)
dim(regulon)



##----------------------------------------------------------------------------##
##                             NMF by cell type                               ##
##----------------------------------------------------------------------------##
library(ArchR)
archrproj <- loadArchRProject(path = params$p_archrproj)
dim(archrproj)


archrproj@sampleMetadata
archrproj@cellColData


unique(samples_df$SampleID)
sort(unique(archrproj@cellColData$Sample))

match(archrproj@cellColData$Sample, samples_df$SampleID_ATAC)
samples_df[!samples_df$SampleID_ATAC %in% archrproj@cellColData$Sample,]

idx <- match(as.character(archrproj@cellColData$Sample), samples_df$SampleID_ATAC)
archrproj@cellColData$months_after_acute_phase <- samples_df$months_after_acute_phase[idx]



as.data.frame(archrproj@cellColData) %>% 
  ggplot(aes(y= predictedGroup)) +
  geom_bar()+
  facet_grid(.~months_after_acute_phase, scales="free") +
  ggtitle("Months after acute phase") +
  cowplot::theme_cowplot()




nmfobj <- readRDS("/home/bq_aquintero/projects/charite_covid19_TF_activity/atac17_rna21/results/integrated/TF_activity_method2_500000/NMF_regulonAUC.RDS")
hm <- HMatrix(nmfobj, k=11)


h <- Heatmap_Regulon(regulon=hm/rowMaxs(hm), seurat=seuratobj, title="H k=11")
dim(hm)
h
