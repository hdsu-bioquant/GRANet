#' Quantification of regulon activity and binarization
#'
#'
#' @param GRANetObject GRANet object with computed regulons.
#' @param min_cell_per_interaction Minimum percentage of cells that are required
#' to have one regulatory interaction. Regulatory interactions found in fewer
#' cells than min_cell_per_interaction are filtered out.
#' @param test_p_val_threshold p-value threshold in order ot keep only regulatory
#' interactions that are related with a change in the target gene expression.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' granetobj <- infer_reg_interactions(GRANetObject=granetobj,
#' min_cell_per_interaction = 0.01,
#' test_p_val_threshold     = 0.05)
#' }
infer_reg_interactions <- function(
  GRANetObject,
  min_cell_per_interaction = 0.01,
  test_p_val_threshold     = 0.05,
  test_by_bins = TRUE,
  binsize      = 1000
){

  #-----------------------------------------------#
  # Make binary matrix of regulatory interactions #
  #-----------------------------------------------#
  message("Making binary matrix of regulatory interactions")
  # Make bin matrix of motif-target x cell
  # PeakID column from motifs granges to subset binary peak counts matrix

  #peak_counts
  motif_target_bin <- Seurat::GetAssayData(
    object = GRANetObject@SeuratObject,
    assay  = GRANetObject@ProjectMetadata$ATAC_assay,
    slot   = "counts")
  motif_target_bin <- as(motif_target_bin, "lgCMatrix")

  # Make bin matrix of motif-TF x cell
  # TF column from motifs granges to subset binary TF activity matrix
  motif_TF_bin <- GRANetObject@RegulonsAUCell$AUCellBin

  # Get expression matrix
  DefaultAssay(GRANetObject@SeuratObject) <- GRANetObject@ProjectMetadata$RNA_assay
  n <- nrow(
    Seurat::GetAssayData(
      object = GRANetObject@SeuratObject,
      assay  = GRANetObject@ProjectMetadata$RNA_assay,
      slot   = "data"))
  if (n == 0) {
    GRANetObject@SeuratObject <- NormalizeData(
      GRANetObject@SeuratObject,
      normalization.method = "LogNormalize",
      scale.factor = 10000)
  }
  exprs_mat <- Seurat::GetAssayData(
    object = GRANetObject@SeuratObject,
    assay  = GRANetObject@ProjectMetadata$RNA_assay,
    slot   = "data")

  # Add regulatory ID
  regID <- GRANetObject@TFmotif_location$annotated_motifs %>%
    as_tibble() %>%
    dplyr::mutate(RegInteractionID = paste0(
      seqnames, "_", start, "_", end, "_",
      TF, "_", target
    )) %>%
    dplyr::select(RegInteractionID)
  GRANetObject@TFmotif_location$annotated_motifs$RegInteractionID <- regID$RegInteractionID


  # Iterate over all bins
  idx <- 1:length(GRANetObject@TFmotif_location$annotated_motifs)
  idx_list <- split(idx, ceiling(seq_along(idx)/binsize))


  regInteractions_list <- pbapply::pblapply(idx_list, function(idx){

    test_reg_interaction(
      motif_target_bin = motif_target_bin, #
      motif_TF_bin     = motif_TF_bin, #
      exprs_mat        = exprs_mat,
      annotated_motifs = GRANetObject@TFmotif_location$annotated_motifs[idx],
      min_cell_per_interaction = 0.01,
      test_p_val_threshold     = 0.05
    )

  })



  reg_inter_bin <- do.call(rbind, lapply(regInteractions_list, "[[", "RegInteractionsBin"))
  reg_rel_test <- do.call(rbind, lapply(regInteractions_list, "[[", "RegInteractionsTest"))


  GRANetObject@Reg_interactions <- list(
    RegInteractionsBin  = reg_inter_bin,
    RegInteractionsTest = reg_rel_test
    #RegInteractionsGRanges = reg_rel_ranges,
  )

  return(GRANetObject)
}





test_reg_interaction <- function(
  motif_target_bin,
  motif_TF_bin,
  exprs_mat,
  annotated_motifs,
  min_cell_per_interaction = 0.01,
  test_p_val_threshold     = 0.05
){

  #-----------------------------------------------#
  # Make binary matrix of regulatory interactions #
  #-----------------------------------------------#
  # message("Making binary matrix of regulatory interactions")
  # Make bin matrix of motif-target x cell
  # PeakID column from motifs granges to subset binary peak counts matrix

  #peak_counts
  motif_target_bin <- motif_target_bin[annotated_motifs$PeakID, , drop=FALSE]

  # Make bin matrix of motif-TF x cell
  # TF column from motifs granges to subset binary TF activity matrix
  motif_TF_bin <- motif_TF_bin[annotated_motifs$RegulonID, , drop=FALSE]

  reg_inter_bin <- motif_target_bin * motif_TF_bin

  #print(reg_inter_bin)

  rownames(reg_inter_bin) <- annotated_motifs$RegInteractionID


  #---------------------------------------#
  # Filtering non-consistent interactions #
  #---------------------------------------#
  # message("Filtering non-consitent interactions")
  idx <- (Matrix::rowSums(reg_inter_bin)/ncol(reg_inter_bin)) >= min_cell_per_interaction
  if (sum(idx) == 0) {
    return(NULL)
  }
  reg_inter_bin <- reg_inter_bin[idx, , drop=FALSE]
  reg_inter_bool <- as(reg_inter_bin, "lgCMatrix")
  #print(rownames(reg_inter_bin))
  #xxxx
  #-----------------------------------------------------#
  # Building matrices to test change in gene expression #
  #-----------------------------------------------------#
  # message("Building matrices to test change in gene expression")
  # Test if the regulatory interaction is related with a change in the
  # expression of he target gene

  # 1. Get annotated motif location of regulatory interactions
  idx <- match(rownames(reg_inter_bin), annotated_motifs$RegInteractionID)
  reg_rel_ranges <- annotated_motifs[idx]

  # 2. Build matrices of gene expression
  exprs_mat <- exprs_mat[reg_rel_ranges$target, , drop=FALSE]
  rownames(exprs_mat) <- rownames(reg_inter_bin)
  exprs_mat <- as.matrix(exprs_mat)

  # 3. Matrix with True when interaction present NA otherwise
  # reg_inter_boolTRUE <- reg_inter_bool
  reg_inter_boolTRUE <- as(reg_inter_bool, "lgeMatrix")
  reg_inter_boolTRUE[reg_inter_boolTRUE == FALSE] <- NA

  # 4. Matrix with True when interaction NOT present NA otherwise
  reg_inter_boolFALSE <- !reg_inter_bool
  reg_inter_boolFALSE[reg_inter_boolFALSE == FALSE] <- NA

  # 5. Multiple indices' matrices with gene expression matrix
  # as matrix is needed to use matrixTests::row_wilcoxon_twosample
  exprs_mat_RegTRUE  <- as.matrix(exprs_mat * reg_inter_boolTRUE)
  exprs_mat_RegFALSE <- as.matrix(exprs_mat * reg_inter_boolFALSE)

  # 6. Split by positive and negative regulatory interactions
  exprs_mat_RegTRUE <- list(
    positive = exprs_mat_RegTRUE[reg_rel_ranges$regulation == 1, , drop=FALSE],
    negative = exprs_mat_RegTRUE[reg_rel_ranges$regulation == -1, , drop=FALSE]
  )
  exprs_mat_RegFALSE <- list(
    positive = exprs_mat_RegFALSE[reg_rel_ranges$regulation == 1, , drop=FALSE],
    negative = exprs_mat_RegFALSE[reg_rel_ranges$regulation == -1, , drop=FALSE]
  )

  #-----------------------------------#
  # Testing change in gene expression #
  #-----------------------------------#
  # message("Testing if the regulatory interaction is related with",
  #         " a change in the target gene expression")

  reg_rel_test <- list(
    positive = matrixTests::row_wilcoxon_twosample(x = exprs_mat_RegTRUE$positive,
                                                   y = exprs_mat_RegFALSE$positive,
                                                   alternative = "greater") %>%
      rownames_to_column("RegInteractionID") %>%
      mutate(regulation = "positive", .after = "RegInteractionID"),
    negative = matrixTests::row_wilcoxon_twosample(x = exprs_mat_RegTRUE$negative,
                                                   y = exprs_mat_RegFALSE$negative,
                                                   alternative = "less") %>%
      rownames_to_column("RegInteractionID") %>%
      mutate(regulation = "negative", .after = "RegInteractionID")

  ) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(pvalue < test_p_val_threshold)


  #---------------------------------------------------------#
  # Retur regulatory interactions to GRANetMultiome object  #
  #---------------------------------------------------------#
  if (nrow(reg_rel_test) == 0) {
    return(NULL)
  }
  idx <- match(reg_rel_test$RegInteractionID, rownames(reg_inter_bin))
  reg_inter_bin <- reg_inter_bin[idx, , drop=FALSE]

  return(list(
    RegInteractionsBin  = reg_inter_bin,
    RegInteractionsTest = reg_rel_test
  ))
}
