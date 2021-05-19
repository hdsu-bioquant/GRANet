# Make sure that TF names are the same used in ArchR
# PeakMatrix
# getPositions
# getPeakSet


#' Title
#'
#' @param GRANetObject
#' @param ArchRProjectObj
#' @param cssClusterArchR
#' @param promoter_size
#' @param cellsWithPeak Minimum percentage of cells in a cluster, that are
#' required to have at least one insertion in a peak to keep it
#' @param threads
#'
#' @return
#' @export
#'
#' @examples
add_motifs_position_from_ArchR <- function(GRANetObject, ArchRProjectObj, cssClusterArchR, promoter_size=5000, cellsWithPeak=0.01, threads=1){

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

  # Get peaks by cell type
  cellTypes <- sort(unique(ArchRProjectObj@cellColData[, cssClusterArchR]))
  names(cellTypes) <- cellTypes
  peaksGr_CellTypes <- lapply(cellTypes, function(cellType){
    # Keep peak if at least one insertion in x% of cells from group
    idx <- BiocGenerics::which(ArchRProjectObj@cellColData[, cssClusterArchR] == cellType)
    peakMatrixCellType <- peakMatrix[,idx]
    percentCells       <- ncol(peakMatrixCellType)*cellsWithPeak
    nCellWithPeak      <- Matrix::rowSums(peakMatrixCellType > 0)
    peaksGr[nCellWithPeak > percentCells]
  })

  message("Computing list of motifs position by cell type...")
  # Create list of motif position by cell type
  motifPositions_Cell <- mclapply(peaksGr_CellTypes, function(peaksGr_Cell){
    endoapply(motifPositions, subsetByOverlaps, peaksGr_Cell)
  }, mc.cores = threads)

  # Add motif position to slot
  GRANetObject@TFmotif_location <- motifPositions_Cell
  GRANetObject@ProjectMetadata$coexprsModulesTFs <- coexprsModulesTFs

  return(GRANetObject)


}
#environment(add_motifs_position_from_ArchR) <- asNamespace('GRANet')

#' Title
#'
#' @param GRANetObject
#' @param promoter_size
#' @param min_regulon_size
#' @param threads
#'
#' @return
#' @export
#'
#' @examples
make_cssRegulons <- function(GRANetObject, promoter_size=5000, min_regulon_size=20, threads=1){

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

  # Remove modules where the TF was not found in the scATAC-seq data
  Coexprs_modules <- GRANetObject@Coexprs_modules %>%
    dplyr::filter(TF %in% GRANetObject@ProjectMetadata$coexprsModulesTFs$TF)

  # Find genes in co-expression modules
  genesGr <- genesGr[genesGr$symbol %in% unique(Coexprs_modules$target)]

  # get TSS
  tssGr <- resize(genesGr, width = 1, fix = "start")

  # Find genes with motif in promoter
  promoters <- trim(suppressWarnings(tssGr+promoter_size))


  #---------------------------------------------------#
  # Filter TFs that have no motif in the ATACseq data #
  #---------------------------------------------------#

  genesMotif_Cell <- lapply(GRANetObject@TFmotif_location, function(motifPositions_Cellgl){
    genesMotif_l <- mclapply(setNames(names(motifPositions_Cellgl), names(motifPositions_Cellgl)), function(TF){
      genesWithMotif <- subsetByOverlaps(promoters, motifPositions_Cellgl[[TF]])
      tibble(TF = TF, target = genesWithMotif$symbol)
    }, mc.cores = threads)

    bind_rows(genesMotif_l) %>%
      dplyr::mutate(TF = sub("_.*", "", TF)) %>%
      dplyr::distinct()
  })

  #-----------------------------------------------#
  # Filter co-expression modules to make regulons #
  #--------------------------------------------- -#
  # Trim out genes that are not a cis target of the TF

  regulons_df_Cell <- lapply(genesMotif_Cell, function(genesMotif){
    inner_join(genesMotif, Coexprs_modules, by=c("TF", "target")) %>%
      dplyr::group_by(TF, regulation) %>%
      dplyr::filter(n() >= min_regulon_size)
  })

  # Add cssRegulons to slot
  GRANetObject@cssRegulons <- regulons_df_Cell

  return(GRANetObject)

}
#environment(make_cssRegulons) <- asNamespace('GRANet')
#make_cssRegulons(GRANetObject = granetobj, promoter_size = 5000, min_regulon_size=20)
#lapply(regulons_df_Cell, function(x) length(unique(x$TF)))


