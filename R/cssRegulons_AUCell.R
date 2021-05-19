
cssRegulons_activity <- function(GRANetObject, threads=1){

  #------------------------------------#
  # Make rankings from gene expression #
  #------------------------------------#

  # Filter cells for which cssRegulon is not available
  idx <- GRANetObject@SeuratObject$cssCluster %in% names(GRANetObject@cssRegulons)


  aucellRankings <- AUCell::AUCell_buildRankings(Seurat::GetAssayData(object = GRANetObject@SeuratObject,
                                                                      slot = "counts")[,idx],
                                                 nCores=threads,
                                                 plotStats=FALSE,
                                                 verbose=TRUE)

  #------------------------------------#
  #           Regulons as list         #
  #------------------------------------#
  # make list of regulons
  regulons_list_Cell <- lapply(GRANetObject@cssRegulons, function(regulons_df){
    regulons_df %>%
      dplyr::mutate(RegulonID = paste0(TF, " (",
                                       if_else(regulation == 1, '+', "-"),
                                       ")")) %>%
      named_group_split(RegulonID) %>%
      # add TF to regulon
      purrr::map(function(x) c(unique(x$TF), x$target))

  })

  #------------------------------------------------------------------------------#
  #               Calculate AUC for each regulon in each cell                    #
  #------------------------------------------------------------------------------#
  # Calculate AUC
  regulonAUC_Cell <- lapply(setNames(names(regulons_list_Cell), names(regulons_list_Cell)), function(CellType){

    CellIDs <- colnames(GRANetObject@SeuratObject)[GRANetObject@SeuratObject$cssCluster %in% CellType]

    AUCell_calcAUC(regulons_list_Cell[[CellType]], aucellRankings[,CellIDs],
                   aucMaxRank = aucellRankings@nGenesDetected["1%"],
                   nCores     = threads)
  })

  regulonAUC_Cell_df <- lapply(regulonAUC_Cell, function(regulonAUC_byCell){
    as.data.frame(regulonAUC_byCell@assays@data$AUC) %>%
      tibble::rownames_to_column("Regulon")
  })

  regulonAUC <- regulonAUC_Cell_df %>%
    purrr::reduce(full_join, by = "Regulon") %>%
    tibble::column_to_rownames("Regulon")


  regulonAUC[is.na(regulonAUC)] <- 0
  regulonAUC <- as.matrix(regulonAUC)
  regulonAUC <- regulonAUC[rowSums(regulonAUC) > 0,]

  return(regulonAUC)





}
#cssRegulons_activity(GRANetObject=granetobj, threads=1)

# rlang
# purrr
# tibble
# Function to split df by groups into a named list
named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = " / ")))

  grouped %>%
    dplyr::group_split() %>%
    rlang::set_names(names)
}
