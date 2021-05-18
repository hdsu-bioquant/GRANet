
# WholeBrainA_62216:    /home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scatac/WholeBrainA_62216_barcoded.bam.arrow
# Lung2_62216:          /home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scatac/Lung2_62216_barcoded.bam.arrow
# Kidney_62016:         /home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scatac/Kidney_62016_barcoded.bam.arrow
# BoneMarrow_62016:     /home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scatac/BoneMarrow_62016_barcoded.bam.arrow
# Liver_62016:          /home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scatac/Liver_62016_barcoded.bam.arrow
# Thymus_62016:         /home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scatac/Thymus_62016_barcoded.bam.arrow
# Spleen_62016:         /home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scatac/Spleen_62016_barcoded.bam.arrow
# Testes_62016:         /home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scatac/Testes_62016_barcoded.bam.arrow
# SmallIntestine_62816: /home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scatac/SmallIntestine_62816_barcoded.bam.arrow

#------------------------------------------------------------------------------#
p_archrproj <- "/home/bq_aquintero/projects/GRANet/mouse_atlas_8tissues/data/scatac/archr/ArchR01_transferLabels/"
#p_seuratobj <- "/home/bq_aquintero/projects/charite_covid19_TF_activity/atac17_rna57/results/scrna/SCENIC/rnaseq_counts.RDS"
annotCol    <- "tissue"
genome      <- "mm9"
Ncores      <- 20
ArrowFiles  <- c("/home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scatac/WholeBrainA_62216_barcoded.bam.arrow",
                 "/home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scatac/Lung2_62216_barcoded.bam.arrow",
                 "/home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scatac/Kidney_62016_barcoded.bam.arrow",
                 "/home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scatac/BoneMarrow_62016_barcoded.bam.arrow",
                 "/home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scatac/Liver_62016_barcoded.bam.arrow",
                 "/home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scatac/Thymus_62016_barcoded.bam.arrow",
                 "/home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scatac/Spleen_62016_barcoded.bam.arrow",
                 "/home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scatac/SmallIntestine_62816_barcoded.bam.arrow")
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

dir.create(dirname(p_archrproj), recursive = TRUE)

# Creat An ArchRProject
archrproj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = p_archrproj,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
archrproj
getAvailableMatrices(archrproj)



#------------------------------------------------------------------------------#
#                      Filter low quality cells                                #
#------------------------------------------------------------------------------#

# log10(min(archrproj$nFrags))
# log10(10000)

archrproj2 <- archrproj[log10(archrproj$nFrags) > 3.0, ]
hist(log10(archrproj$nFrags))
hist(log10(archrproj2$nFrags))

archrproj
archrproj2

archrproj2 <- my_filterDoublets(archrproj2)
hist(log(archrproj2$nFrags))
archrproj2


# 4.2 Iterative Latent Semantic Indexing (LSI)
archrproj2 <- addIterativeLSI(
  ArchRProj = archrproj2,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(5),
    sampleCells = 10000,
    maxClusters = 20,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30,
  UMAPParams = list(n_neighbors = 120, min_dist = 0.25)
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


cM <- confusionMatrix(paste0(archrproj2$Clusters), paste0(archrproj2$Sample))
cM

cM <- as.matrix(cM / Matrix::rowSums(cM))
p <- Heatmap(matrix = as.matrix(cM),
             col = viridis(100))
pdf(file.path(p_archrproj, "/clusters_heatmap.pdf"))
p
dev.off()


#------------------------------------------------------------------------------#
#                                6 UMAP                                        #
#------------------------------------------------------------------------------#
# 6.1 Uniform Manifold Approximation and Projection (UMAP)
archrproj2 <- addUMAP(
  ArchRProj = archrproj2,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 80,
  minDist = 0.25,
  metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = archrproj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = archrproj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
#p1 + p2
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = archrproj2,
        addDOC = FALSE, width = 5, height = 5)

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
#                              Add tissue identity                             #
#------------------------------------------------------------------------------#
head(archrproj2@cellColData)
archrproj2$tissue <- paste0("Adult", gsub("_.*|Whole|2|A", "", archrproj2$Sample))

p1 <- plotEmbedding(archrproj2, colorBy = "cellColData", name = "tissue")

plotPDF(p1, name = "Plot-UMAP-Remap-Clusters.pdf", ArchRProj = archrproj2, addDOC = FALSE, width = 5, height = 5)

h <- Heatmap(as.matrix(cM), col=magma(100), name="Confusion Matrix\nCluster vs Tissue")
plotPDF(h, name = "Confusion-Matrix-Cluster-vs-Tissue.pdf", ArchRProj = archrproj2, addDOC = FALSE, width = 7, height = 7)

saveArchRProject(ArchRProj = archrproj2, load = FALSE)


#------------------------------------------------------------------------------#
#                   Chapter 9 Pseudo-bulk Replicates in ArchR                  #
#------------------------------------------------------------------------------#
archrproj4 <- addGroupCoverages(ArchRProj = archrproj2, groupBy = "tissue")


#------------------------------------------------------------------------------#
#                   Chapter 10 Calling Peaks with ArchR                        #
#------------------------------------------------------------------------------#
pathToMacs2 <- findMacs2()
#pathToMacs2 <- "/home/bq_aquintero/miniconda3/archerhelper/bin/macs2"

archrproj4 <- addReproduciblePeakSet(
  ArchRProj = archrproj4,
  groupBy = "tissue",
  pathToMacs2 = pathToMacs2
)


archrproj5 <- addPeakMatrix(archrproj4)
getAvailableMatrices(archrproj5)
ArchR::getPeakSet(archrproj5)

#------------------------------------------------------------------------------#
#            Chapter 12 Motif and Feature Enrichment with ArchR                #
#------------------------------------------------------------------------------#
# 12.1 Motif Enrichment in Differential Peaks
archrproj5 <- addMotifAnnotations(ArchRProj = archrproj5, motifSet = "cisbp", name = "Motif", force=TRUE)
getPositions(archrproj5)



#saveArchRProject(ArchRProj = archrproj5, load = FALSE)
saveArchRProject(ArchRProj = archrproj5,
                 outputDirectory = "/home/bq_aquintero/projects/GRANet/mouse_atlas_8tissues/data/scatac/archr/ArchR02_MotifMatch/",
                 load = FALSE)

