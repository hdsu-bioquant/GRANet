options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

p_archrproj    <- dirname(as.character(args[1]))
p_archrprojout <- dirname(as.character(args[2]))
Ncores         <- as.numeric(args[3])
#------------------------------------------------------------------------------#
# p_archrproj    <- "~/projects/mouse_atlas_TF_activity/adult_9tissues/results/scatac/archr/ArchR01_transferLabels/"
# p_archrprojout <- "~/projects/mouse_atlas_TF_activity/adult_9tissues/results/scatac/archr/ArchR02_MotifMatch/"
# Ncores      <- 20
#------------------------------------------------------------------------------#
p_archrproj    <- "/home/bq_aquintero/projects/charite_covid19_TF_activity/atac17_rna21/results/scatac/archr/ArchR01_transferLabels"
p_archrprojout <- "/home/bq_aquintero/projects/charite_covid19_TF_activity/atac17_rna21/results/scatac/archr/ArchR02_MotifMatch/"
Ncores      <- 20

library(ArchR)
library(patchwork)
library(ComplexHeatmap)
library(viridis)

ArchR::addArchRThreads(Ncores)
#addArchRGenome(genome)


#------------------------------------------------------------------------------#
#                          Debugged functions                                  #
#------------------------------------------------------------------------------#

archrproj3 <- loadArchRProject(path = p_archrproj)

archrproj3$Clusters2 <- archrproj3$predictedGroup

dim(archrproj3)
getGenome(archrproj3)
#saveArchRProject(ArchRProj = archrproj3, outputDirectory = p_archrprojout, load = FALSE)

#archrproj3@peakAnnotation
#------------------------------------------------------------------------------#
#                   Chapter 9 Pseudo-bulk Replicates in ArchR                  #
#------------------------------------------------------------------------------#
archrproj4 <- addGroupCoverages(ArchRProj = archrproj3, groupBy = "Clusters2")


#------------------------------------------------------------------------------#
#                   Chapter 10 Calling Peaks with ArchR                        #
#------------------------------------------------------------------------------#
pathToMacs2 <- findMacs2()
#pathToMacs2 <- "/home/bq_aquintero/miniconda3/archerhelper/bin/macs2"

archrproj4 <- addReproduciblePeakSet(
  ArchRProj = archrproj4, 
  groupBy = "Clusters2", 
  pathToMacs2 = pathToMacs2
)

getPeakSet(archrproj4)

#saveArchRProject(ArchRProj = archrproj4, outputDirectory = "Save-archrproj4", load = FALSE)
archrproj5 <- addPeakMatrix(archrproj4)
getAvailableMatrices(archrproj5)
ArchR::getPeakSet(archrproj5)
# #------------------------------------------------------------------------------#
# #             Chapter 11 Identifying Marker Peaks with ArchR                   #
# #------------------------------------------------------------------------------#
# #11.1 Identifying Marker Peaks with ArchR
# table(archrproj5$Clusters2)
# markersPeaks <- getMarkerFeatures(
#   ArchRProj = archrproj5, 
#   useMatrix = "PeakMatrix", 
#   groupBy = "Clusters2",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# markersPeaks
# 
# markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
# markerList
# 
# 
# heatmapPeaks <- markerHeatmap(
#   seMarker = markersPeaks, 
#   cutOff = "FDR <= 0.05 & Log2FC >= 3",
#   transpose = TRUE
# )
# 
# dim(heatmapPeaks@matrix)
# 
# plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = archrproj5, addDOC = FALSE)
# 

#------------------------------------------------------------------------------#
#            Chapter 12 Motif and Feature Enrichment with ArchR                #
#------------------------------------------------------------------------------#
# 12.1 Motif Enrichment in Differential Peaks
#archrproj5 <- loadArchRProject(path = p_archrproj)
archrproj5 <- addMotifAnnotations(ArchRProj = archrproj5, motifSet = "cisbp", name = "Motif", force=TRUE)
#archrproj5@peakAnnotation
# # 12.3 ArchR Enrichment
# # 12.3.1 Encode TF Binding Sites
# archrproj5 <- addArchRAnnotations(ArchRProj = archrproj5, collection = "EncodeTFBS")
# enrichEncode <- peakAnnoEnrichment(
#   seMarker = markersPeaks,
#   ArchRProj = archrproj5,
#   peakAnnotation = "EncodeTFBS",
#   cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
# )
# heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)
# plotPDF(heatmapEncode, name = "10_EncodeTFBS-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = archrproj5, addDOC = FALSE)


# archrproj5 <- my_addArchRAnnotations(ArchRProj = archrproj5, collection = "ATAC")
# enrichATAC <- peakAnnoEnrichment(
#   seMarker = markersPeaks,
#   ArchRProj = archrproj5,
#   peakAnnotation = "ATAC",
#   cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
# )
# heatmapATAC <- plotEnrichHeatmap(enrichATAC, n = 7, transpose = TRUE)
# plotPDF(heatmapATAC, name = "11_ATAC-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = archrproj5, addDOC = FALSE)


# archrproj5 <- my_addArchRAnnotations(ArchRProj = archrproj5, collection = "Codex")
# enrichCodex <- peakAnnoEnrichment(
#   seMarker = markersPeaks,
#   ArchRProj = archrproj5,
#   peakAnnotation = "Codex",
#   cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
# )
# heatmapCodex <- plotEnrichHeatmap(enrichCodex, n = 7, transpose = TRUE)
# plotPDF(heatmapCodex, name = "12_Codex-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = archrproj5, addDOC = FALSE)
# 
getPositions(archrproj5)
#------------------------------------------------------------------------------#
#            Chapter 13 ChromVAR Deviatons Enrichment with ArchR               #
#------------------------------------------------------------------------------#
# 13.1 Motif Deviations
archrproj5 <- addBgdPeaks(archrproj5)
archrproj5 <- addDeviationsMatrix(
  ArchRProj = archrproj5, 
  #peakAnnotation = "Motif",
  force = TRUE
)
#archrproj5@peakAnnotation
plotVarDev <- getVarDeviations(archrproj5, name = "MotifMatrix", plot = TRUE)
#plotVarDev
plotPDF(plotVarDev, name = "13_Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = archrproj5, addDOC = FALSE)


#saveArchRProject(ArchRProj = archrproj5, load = FALSE)
saveArchRProject(ArchRProj = archrproj5, outputDirectory = p_archrprojout, load = FALSE)

