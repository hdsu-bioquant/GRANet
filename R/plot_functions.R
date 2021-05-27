
#' Heatmap of cssRegulons activities
#'
#' @param GRANetObject GRANet object with calculated cssRegulon activities.
#' See function "cssRegulons_activity".
#' @param show_column_names Display cell IDs.
#' @param show_row_names  Display transcription factor symbol.
#' @param show_column_dend Display column dendrogram.
#' @param show_row_dend Display row dendrogram
#' @param cluster_rows Cluster rows (cssRegulons).
#'
#' @return
#' @export
#'
#' @examples
Heatmap_AUCell <- function(GRANetObject,
                           show_column_names = FALSE,
                           show_row_names    = FALSE,
                           show_column_dend  = TRUE,
                           show_row_dend     = TRUE,
                           cluster_rows      = TRUE){
  heat_anno <- ComplexHeatmap::HeatmapAnnotation(
    df= data.frame(cssCluster = GRANetObject@SeuratObject$cssCluster))

  ComplexHeatmap::Heatmap(GRANetObject@cssRegulonsAUCell,
          col = viridis::viridis(100),
          name = "cssRegulon\nAUCell",
          top_annotation    = heat_anno,
          show_column_names = show_column_names,
          show_row_names    = show_row_names,
          show_column_dend  = show_column_dend,
          show_row_dend     = show_row_dend,
          cluster_rows      = cluster_rows)

}
