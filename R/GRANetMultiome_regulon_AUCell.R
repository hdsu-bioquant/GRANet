#' Quantification of regulon activity and binarization
#'
#'
#' @param GRANetObject GRANet object with computed regulons.
#' @param aucMaxRank possible values: "min", "1\%", "5\%", "10\%", "50\%", "100\%".
#' @param threads Number of threads to use.
#' @param fastBinarization Regulon activity binarization can be run in a
#' multithreaded implementation. However, for large dataset the single threaded
#' implementation is more stable.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' granetobj <- regulons_activity(GRANetObject=granetobj,
#' aucMaxRank="50%", threads=1)
#' }
regulons_activity <- function(
  GRANetObject,
  aucMaxRank       = "1%",
  threads          = 1,
  fastBinarization = FALSE
  ){

  message("Building rankind and computing TF activity")
  #------------------------------------#
  # Make rankings from gene expression #
  #------------------------------------#
  aucellRankings <- AUCell::AUCell_buildRankings(
    exprMat   = Seurat::GetAssayData(
      object = GRANetObject@SeuratObject,
      assay  = GRANetObject@ProjectMetadata$RNA_assay,
      slot   = "counts"),
    nCores    = threads,
    plotStats = FALSE,
    verbose   = TRUE)

  #---------------------------------------------#
  # Calculate AUC for each regulon in each cell #
  #---------------------------------------------#
  regulonAUC_Cell <- AUCell::AUCell_calcAUC(
    geneSets   = GRANetObject@Regulons,
    rankings   = aucellRankings,
    aucMaxRank = aucellRankings@nGenesDetected[aucMaxRank],
    nCores     = threads)

  regulonAUC <- regulonAUC_Cell@assays@data$AUC


  regulonAUC[is.na(regulonAUC)] <- 0
  regulonAUC <- as.matrix(regulonAUC)
  regulonAUC <- regulonAUC[rowSums(regulonAUC) > 0,]
  regulonAUC <- as.data.frame(t(regulonAUC))


  #-----------------------------------------------#
  # Binarize regulon AUC into active/inactive TFs #
  #-----------------------------------------------#
  message("Binarizing TF activity...")
  pyscenic <- reticulate::import("pyscenic")

  if (fastBinarization) {
    scenic_binarized <- binarize_AUCellScores_fast(regulonAUC)
  } else {
    scenic_binarized <- binarize_AUCellScores_slow(regulonAUC)
  }

  # Add cssRegulons AUCell to slot
  GRANetObject@RegulonsAUCell <- list(
    AUCellScores = regulonAUC,
    AUCellBin    = scenic_binarized$AUCellBin,
    AUCellThrs   = scenic_binarized$AUCellThrs
  )


  return(GRANetObject)

}

binarize_AUCellScores_fast <- function(AUCellScores){
  pyscenic <- reticulate::import("pyscenic")

  pyscenic <- reticulate::import("pyscenic")
  scenic_binarized <- pyscenic$binarization$binarize(
    auc_mtx     = AUCellScores,
    num_workers = as.integer(threads))

  return(list(AUCellBin  = as(t(scenic_binarized[[1]]), "lgCMatrix"),
              AUCellThrs = scenic_binarized[[2]]))
}


binarize_AUCellScores_slow <- function(AUCellScores){
  pyscenic <- reticulate::import("pyscenic")

  thresholds <- pbapply::pblapply(colnames(AUCellScores), function(regulonID){
    pyscenic$binarization$derive_threshold(auc_mtx      = AUCellScores,
                                           regulon_name = regulonID )
  })
  thresholds <- unlist(thresholds)
  AUCellBin <- (t(AUCellScores) > thresholds) * 1
  AUCellBin <- as(AUCellBin, "lgCMatrix")

  return(list(AUCellBin  = AUCellBin,
              AUCellThrs = thresholds))
}

