# Make sure that TF names are the same used in ArchR
granetobj
granetobj@ProjectMetadata$Genome <- "mm9"

PeakMatrix
getPositions
getPeakSet

library(dplyr)

archrproj@cellColData[, "tissue"]

# cellWithPeak minimum percentage of cells in a cluster, that are required to have at least one insertion in a peak to keep it

extract_motifs_ArchR <- function(GRANetObject, ArchRProjectObj, cssClusterArchR, promoter_size=5000, cellsWithPeak=0.01, threads=1){
  genome <- GRANetObject@ProjectMetadata$Genome

  #----------------------------------------------#
  # Find genes part of the co-expression modules #
  #----------------------------------------------#
  message("Extracting genomic location genes in co-expression modules...")
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
  message("Filtering TFs from co-expression modules not found in the ArchR project peak annotation position (motifs)...")
  motifPositions <- ArchR::getPositions(ArchRProjectObj)

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


  #-------------------------------------------#
  # Find genes with a motif in their promoter #
  #-------------------------------------------#
  message("Extracting peak matrix from ArchR project...")
  # Find peaks around the TSS of the target genes
  # To identify genes that are cis-regulated by the TFs
  # And later remove genes that are downstream target of the TF

  # Split peaks into peaks found by Cell type

  # Load peaks for each cell type
  peaksGr <- ArchR::getPeakSet(ArchRProjectObj)

  # Get peaks matrix
  peakMatrix <- ArchR::getMatrixFromProject(ArchRProjectObj, useMatrix = "PeakMatrix", logFile=NULL)
  peakMatrix <- SummarizedExperiment::assay(peakMatrix)

  #return(peakMatrix)

  # Get peaks by cell type
  cellTypes <- sort(unique(ArchRProjectObj@cellColData[, cssClusterArchR]))
  names(cellTypes) <- cellTypes
  peaksGr_CellTypes <- lapply(cellTypes, function(cellType){
    # Keep peak if at least one insertion in x% of cells from group
    idx <- which(ArchRProjectObj@cellColData[, cssClusterArchR] == cellType)
    peakMatrixCellType <- peakMatrix[,idx]
    percentCells       <- ncol(peakMatrixCellType)*cellsWithPeak
    nCellWithPeak      <- rowSums(peakMatrixCellType > 0)
    peaksGr[nCellWithPeak > percentCells]
  })

  peaksGr_CellTypes

  message("Computing list of motifs position by cell type...")
  # Create list of motif position by cell type
  motifPositions_Cell <- mclapply(peaksGr_CellTypes, function(peaksGr_Cell){
    endoapply(motifPositions, subsetByOverlaps, peaksGr_Cell)
  }, mc.cores = threads)

  # Find genes with motif in promoter
  promoters <- trim(tssGr+promoter_size)

  genesMotif_Cell <- lapply(motifPositions_Cell, function(motifPositions_Cellgl){
    genesMotif_l <- mclapply(setNames(names(motifPositions_Cellgl), names(motifPositions_Cellgl)), function(TF){
      genesWithMotif <- subsetByOverlaps(promoters, motifPositions_Cellgl[[TF]])
      tibble(TF = TF, target = genesWithMotif$symbol)
    }, mc.cores = threads)

    #print(genesMotif_l)
    bind_rows(genesMotif_l) %>%
      dplyr::mutate(TF = sub("_.*", "", TF)) %>%
      dplyr::distinct()
  })

  return(genesMotif_Cell)


}


x <- extract_motifs_ArchR(GRANetObject=granetobj, ArchRProjectObj=archrproj, cssClusterArchR="tissue", threads=8)



lapply(sort(unique(archrproj@cellColData[, "tissue"])), function(cellType){
  # Keep peak if at least one insertion in x% of cells from group
  idx <- which(archrproj@cellColData[, "tissue"] == cellType)
  peakMatrixCellType <- x[,idx]
  percentCells       <- ncol(peakMatrixCellType)*cellsWithPeak
  nCellWithPeak      <- rowSums(peakMatrixCellType > 0)
  peaksGr[nCellWithPeak > percentCells]
})



x[1:5,1:5]

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

