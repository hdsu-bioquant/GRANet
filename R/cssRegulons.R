# Make sure that TF names are the same used in ArchR
granetobj
granetobj@ProjectMetadata$Genome <- "mm9"

library(dplyr)

extract_motifs_ArchR <- function(GRANetObject, ArchRProjectObj){
  genome <- GRANetObject@ProjectMetadata$Genome

  #----------------------------------------------#
  # Find genes part of the co-expression modules #
  #----------------------------------------------#
  if (genome == "hg19") {
    library(EnsDb.Hsapiens.v75)
    genesGr <- genes(EnsDb.Hsapiens.v75)
  } else if (genome == "hg38") {
    library(EnsDb.Hsapiens.v86)
    genesGr <- genes(EnsDb.Hsapiens.v86)
  } else if (genome == "mm9") {
    library(EnsDb.Mmusculus.v75)
    genesGr <- genes(EnsDb.Mmusculus.v75)
  } else if (genome == "mm10") {
    library(EnsDb.Mmusculus.v79)
    genesGr <- genes(EnsDb.Mmusculus.v79)
  }
  seqlevelsStyle(genesGr) <- 'UCSC'

  # Find genes in co-expression modules
  genesGr <- genesGr[genesGr$symbol %in% unique(GRANetObject@Coexprs_modules$target)]

  # get TSS
  tssGr <- resize(genesGr, width = 1, fix = "start")
  #seqlevelsStyle(tssGr) <- 'UCSC'
  tssGr


  #---------------------------------------------------#
  # Filter TFs that have no motif in the ATACseq data #
  #---------------------------------------------------#
  motifPositions <- getPositions(ArchRProjectObj)

  coexprsModulesTFs <- tibble(motif_archr = names(motifPositions),
                      TF_archr    = sub("_.*", "", names(motifPositions))) %>%
    dplyr::mutate(TF = GRANetObject@Coexprs_modules$TF[match(TF_archr, GRANetObject@Coexprs_modules$TF)]) %>%
    dplyr::filter(!is.na(TF)) %>%
    dplyr::distinct()


  # Remove modules where the TF was not found in the scATAC-seq data
  Coexprs_modules <- GRANetObject@Coexprs_modules %>%
    dplyr::filter(TF %in% coexprsModulesTFs$TF)

  # Remove motifs for TF not found
  motifPositions <- motifPositions[names(motifPositions) %in% coexprsModulesTFs$motif_archr]
  length(motifPositions)
  motifPositions
}


extract_motifs_ArchR(GRANetObject=granetobj, ArchRProjectObj=archrproj)


#------------------------------------------------------------------------------#
#                 Find genes with a motif in its promoter.                     #
#------------------------------------------------------------------------------#
# Find peaks between 5k of tss
# To identify genes that are cis regulated by the TF
# And later remove genes that are downstream target of the TF

# Split peaks into peaks found by Cell type
# This way cell type-specific regulos will be buils


# Load peaks for each cell type
peaksGr <- ArchR::getPeakSet(archrproj)
unique(names(peaksGr))

# Get peaks matrix
ArchR::getAvailableMatrices(archrproj)
peakMatrix <- ArchR::getMatrixFromProject(archrproj, useMatrix = "PeakMatrix")
peakMatrix <- assay(peakMatrix)
dim(peakMatrix)
length(peaksGr)

# Get peaks by cell type
#cellWithPeak <- 0.01
table(archrproj$predictedGroup)
cellTypes <- sort(unique(archrproj$predictedGroup))
names(cellTypes) <- cellTypes
peaksGr_CellTypes <- lapply(cellTypes, function(cellType){
  # Keep peak if at least one insertion in x% of cells from group
  peakMatrixCellType <- peakMatrix[,archrproj$predictedGroup == cellType]
  percentCells         <- ncol(peakMatrixCellType)*cellWithPeak
  nCellWithPeak      <- rowSums(peakMatrixCellType > 0)
  peaksGr[nCellWithPeak > percentCells]
})

# Create list of motif position by cell type
motifPositions_Cell <- lapply(peaksGr_CellTypes, function(peaksGr_Cell){
  endoapply(motifPositions, subsetByOverlaps, peaksGr_Cell)
})
#motifPositions_Cell
sapply(motifPositions_Cell, function(x){
  summary(sapply(x, length))
})



# Find genes with motif in promoter
promoters <- trim(tssGr+promoter_size)
#motifPositions_Cell$Ciliated

genesMotif_Cell <- lapply(motifPositions_Cell, function(motifPositions_Cellgl){
  genesMotif_l <- mclapply(setNames(names(motifPositions_Cellgl), names(motifPositions_Cellgl)), function(TF){
    genesWithMotif <- subsetByOverlaps(promoters, motifPositions_Cellgl[[TF]])
    tibble(TF = TF, target = genesWithMotif$symbol)
    #}, mc.cores = 1)
  }, mc.cores = Ncores)

  #print(genesMotif_l)
  bind_rows(genesMotif_l) %>%
    dplyr::mutate(TF = sub("_.*", "", TF)) %>%
    dplyr::distinct()
})

genesMotif_Cell
sapply(genesMotif_Cell, nrow)
#------------------------------------------------------------------------------#
#             Filter co-expression modules to make regulons                    #
#------------------------------------------------------------------------------#
# Trim out genes that are not a cis target of the TF

regulons_df_Cell <- lapply(genesMotif_Cell, function(genesMotif){
  inner_join(genesMotif, scenic_corrmodules, by=c("TF", "target")) %>%
    dplyr::group_by(TF, regulation) %>%
    dplyr::filter(n() >= min_regulon_size)
})
regulons_df_Cell
lapply(regulons_df_Cell, function(x) length(unique(x$TF)))

saveRDS(regulons_df_Cell, file = p_regulons)

