options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

p_archrproj <- dirname(as.character(args[1]))
p_seuratobj <- as.character(args[2])
annotCol    <- as.character(args[3])
genome      <- as.character(args[4])
Ncores      <- as.numeric(args[5])
ArrowFiles  <- as.character(args[6:length(args)])
#------------------------------------------------------------------------------#
# p_archrproj <- "~/projects/mouse_atlas_TF_activity/adult_lung/results/scatac/archr/ArchR01_transferLabels/"
# p_seuratobj <-  "~/projects/mouse_atlas_TF_activity/data/scrna/Han_atlas/MCA_AdultLung_Seurat.RDS"
# annotCol    <- "Annotation"
# genome      <- "mm9"
# Ncores      <- 20
# ArrowFiles  <- c("~/projects/mouse_atlas_TF_activity/data/scatac/Lung2_62216_barcoded.bam.arrow")
# #------------------------------------------------------------------------------#
# p_archrproj <- "/home/bq_aquintero/projects/charite_covid19_TF_activity/atac17_rna21/results/scatac/archr/ArchR01_transferLabels/" 
# p_seuratobj <-  "/home/bq_aquintero/projects/charite_covid19_TF_activity/results/scrna/seurat/rna_seurat_int_annot_transfer.RDS"
# annotCol    <- "BioClassification"
# genome      <- "hg19"
# Ncores      <- 20
# ArrowFiles  <- c("/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/SC2-10025AT-NS02FU.arrow", 
#                  "/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/SC2-30026AT-NS02FU.arrow", 
#                  "/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/SC2-21001AT-NS01FU.arrow", 
#                  "/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/SC2-20025AT-NS02FU.arrow",
#                  "/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/SC2-30011AT-NS02FU.arrow", 
#                  "/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/SC2-21004AT-NS01FU.arrow", 
#                  "/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/SC2-21018AT-NS01FU.arrow", 
#                  "/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/SC2-21013AT-NS01FU.arrow", 
#                  "/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/SC2-21010AT-NS01FU.arrow", 
#                  "/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/SC2-20017AT-NS02FU.arrow", 
#                  "/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/SC2-10021AT-NS02FU.arrow", 
#                  "/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/SC2-10012AT-NS02FU.arrow", 
#                  "/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/SC2-21011AT-NS01FU.arrow", 
#                  "/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/SC2-21016AT-NS01FU.arrow", 
#                  "/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/SC2-10021AT-NS03FU2.arrow", 
#                  "/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/SC2-50069AT-NS01.arrow", 
#                  "/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/SC2-50049AT-NS01.arrow")
#------------------------------------------------------------------------------#

library(ArchR)
library(patchwork)
library(ComplexHeatmap)
library(viridis)

#------------------------------------------------------------------------------#
#                          Debugged functions                                  #
#------------------------------------------------------------------------------#

my_filterDoublets <- function(ArchRProj = NULL, cutEnrich = 1, cutScore = -Inf, filterRatio = 1){
  
  fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
  for (i in seq_along(fn)) {
    tryCatch({
      eval(parse(text = paste0(fn[i], "<-ArchR:::", fn[i])))
    }, error = function(x) {
    })
  }
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = cutEnrich, name = "cutEnrich", valid = c("numeric"))
  .validInput(input = cutScore, name = "cutScore", valid = c("numeric"))
  .validInput(input = filterRatio, name = "filterRatio", valid = c("numeric"))
  
  if(any(grepl("filterDoublets", names(ArchRProj@projectSummary)))){
    stop("Already ran filterDoublets on ArchRProject! Cannot be re-ran on an ArchRProject!")
  }
  
  df <- getCellColData(ArchRProj, c("Sample", "DoubletEnrichment", "DoubletScore"))
  splitDF <- split(seq_len(nrow(df)), as.character(df$Sample))
  
  cellsFilter <- lapply(splitDF, function(y){
    
    x <- df[y, ,drop = FALSE]
    
    n <- nrow(x)
    
    x <- x[order(x$DoubletEnrichment, decreasing = TRUE), ]
    
    if(!is.null(cutEnrich)){
      x <- x[which(x$DoubletEnrichment >= cutEnrich), ]
    } 
    
    if(!is.null(cutScore)){
      x <- x[which(x$DoubletScore >= cutScore), ]
    } 
    
    if(nrow(x) > 0){
      head(rownames(x), filterRatio * n * (n / 100000))
    }else{
      NULL
    }
    
  }) %>% unlist(use.names=FALSE)
  
  message("Filtering ", length(cellsFilter), " cells from ArchRProject!")
  tabRemove <- table(df[cellsFilter,]$Sample)
  tabAll <- table(df$Sample)
  samples <- unique(df$Sample)
  for(i in seq_along(samples)){
    if(!is.na(tabRemove[samples[i]])){
      message("\t", samples[i], " : ", tabRemove[samples[i]], " of ", tabAll[samples[i]], " (", round(100 * tabRemove[samples[i]] / tabAll[samples[i]], 1),"%)")
    }else{
      message("\t", samples[i], " : ", 0, " of ", tabAll[samples[i]], " (0%)")
    }
  }
  
  if(length(cellsFilter) > 0){
    
    ArchRProj@cellColData <- ArchRProj@cellColData[rownames(ArchRProj@cellColData) %ni% cellsFilter,,drop=FALSE]
    
  }
  
  ArchRProj
  
}

#------------------------------------------------------------------------------#
#                                Create ArchR                                  #
#------------------------------------------------------------------------------#




ArchR::addArchRThreads(Ncores)
addArchRGenome(genome)


# Creat An ArchRProject
archrproj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = p_archrproj,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
archrproj
getAvailableMatrices(archrproj)

#saveArchRProject(ArchRProj = archrproj, outputDirectory = "Save-archrproj", load = FALSE)
#archrproj <- loadArchRProject("/media/ag-cherrmann/projects/10_charite_covid19/results/scatac/archr/Save-archrproj/")

archrproj2 <- my_filterDoublets(archrproj)

# 4.2 Iterative Latent Semantic Indexing (LSI)
archrproj2 <- addIterativeLSI(
  ArchRProj = archrproj2,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

#4.4 Batch Effect Correction wtih Harmony
archrproj2 <- addHarmony(
  ArchRProj = archrproj2,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)


# 5.1 Clustering using Seurat’s FindClusters() function
archrproj2 <- addClusters(
  input = archrproj2,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

head(archrproj2$Clusters)
table(archrproj2$Clusters)

if (length(unique(archrproj2$Sample)) > 1) {
  cM <- confusionMatrix(paste0(archrproj2$Clusters), paste0(archrproj2$Sample))
  cM
  
  cM <- as.matrix(cM / Matrix::rowSums(cM))
  p <- Heatmap(matrix = as.matrix(cM), 
               col = viridis(100))
  pdf(file.path(p_archrproj, "/clusters_heatmap.pdf"))
  p
  dev.off()
  
}

#saveArchRProject(ArchRProj = archrproj2, outputDirectory = "Save-archrproj2", load = FALSE)

#------------------------------------------------------------------------------#
#                                6 UMAP                                        #
#------------------------------------------------------------------------------#
# 6.1 Uniform Manifold Approximation and Projection (UMAP)
archrproj2 <- addUMAP(
  ArchRProj = archrproj2, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = archrproj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = archrproj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
#p1 + p2
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = archrproj2, 
        addDOC = FALSE, width = 5, height = 5)

#6.3 Dimensionality Reduction After Harmony
archrproj2 <- addUMAP(
  ArchRProj = archrproj2, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
p3 <- plotEmbedding(ArchRProj = archrproj2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = archrproj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
#p3 + p4
plotPDF(p1,p2,p3,p4, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = archrproj2, addDOC = FALSE, width = 5, height = 5)

#------------------------------------------------------------------------------#
#                            7.3 Identifying Marker Genes                      #
#------------------------------------------------------------------------------#
markersGS <- getMarkerFeatures(
  ArchRProj = archrproj2, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
lapply(markerList, dim)

markerGenes <- do.call(c, lapply(markerList, function(x) x$name[1:2]))

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
#ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = archrproj2, addDOC = FALSE)


#------------------------------------------------------------------------------#
#                8 Defining Cluster Identity with scRNA-seq                    #
#------------------------------------------------------------------------------#
seuratobj <- readRDS(p_seuratobj)
#seuratobj <- NormalizeData(object = seuratobj, verbose = FALSE)
head(seuratobj@meta.data)
DefaultAssay(seuratobj) <- "RNA"

#table(seuratobj$BioClassification)

archrproj2 <- addGeneIntegrationMatrix(
  ArchRProj = archrproj2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seuratobj,
  addToArrow = FALSE,
  groupRNA = annotCol,
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

table(archrproj2$predictedGroup_Un)

p1 <- plotEmbedding(
  archrproj2, 
  colorBy = "cellColData", 
  name = "predictedGroup_Un"
)
#p1
plotPDF(p1, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = archrproj2, addDOC = FALSE, width = 5, height = 5)

#saveArchRProject(ArchRProj = archrproj2, outputDirectory = "Save-archrproj2", load = FALSE)

#8.2 Adding Pseudo-scRNA-seq profiles for each scATAC-seq cell
archrproj3 <- addGeneIntegrationMatrix(
  ArchRProj = archrproj2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seuratobj,
  addToArrow = TRUE,
  force= TRUE,
  #groupList = groupList,
  groupRNA = annotCol,
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)
getAvailableMatrices(archrproj3)
archrproj3 <- addImputeWeights(archrproj3)


#8.3 Labeling scATAC-seq clusters with scRNA-seq information
cM <- confusionMatrix(archrproj3$Clusters, archrproj3$predictedGroup)
labelOld <- rownames(cM)
labelOld
# labelNew <- colnames(cM)[apply(cM, 1, which.max)]
# labelNew

archrproj3$Clusters2 <- archrproj3$predictedGroup

p1 <- plotEmbedding(archrproj3, colorBy = "cellColData", name = "Clusters2")
#p1

plotPDF(p1, name = "Plot-UMAP-Remap-Clusters.pdf", ArchRProj = archrproj3, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = archrproj3, load = FALSE)
