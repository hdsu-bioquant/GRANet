#' Quantification of regulon activity and binarization
#'
#'
#' @param GRANetObject GRANet object with computed regulons.
#' @param aucMaxRank possible values: "min", "1\%", "5\%", "10\%", "50\%", "100\%".
#' @param threads Number of threads to use.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' granetobj <- regulons_activity(GRANetObject=granetobj,
#' aucMaxRank="50%", threads=1)
#' }
regulons_activity <- function(GRANetObject, aucMaxRank="1%", threads=1){

  message("Building rankind and computing TF acivty...")
  #------------------------------------#
  # Make rankings from gene expression #
  #------------------------------------#
  aucellRankings <- AUCell::AUCell_buildRankings(
    exprMat   = Seurat::GetAssayData(object = GRANetObject@SeuratObject, slot = "counts"),
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
  scenic_binarized <- pyscenic$binarization$binarize(
    auc_mtx     = regulonAUC,
    num_workers = as.integer(threads))

  # Add cssRegulons AUCell to slot
  GRANetObject@RegulonsAUCell <- list(
    AUCellScores = regulonAUC,
    AUCellBin    = as(t(scenic_binarized[[1]]), "lgCMatrix"),
    AUCellThrs   = scenic_binarized[[2]]
  )


  return(GRANetObject)

}
