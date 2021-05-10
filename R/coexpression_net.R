


#' Title
#'
#' @slot GeneExpression list.
#' @slot ArchR list.
#' @slot TFmotif_location list.
#' @slot Coexprs_modules data.frame.
#' @slot cssRegulons list.
#' @slot SignFeatures data.frame.
#'
#' @return
#' @export
#'
#' @examples
ButchR_NMF <- setClass(
  Class = "GRANet",
  slots = list(GeneExpression   = "list",
               ArchR            = "list",
               TFmotif_location = "list",
               Coexprs_modules  = "data.frame",
               cssRegulons      = "list",
               SignFeatures = "data.frame" )

)





#' Title
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
gene_coex_net <- function(obj){
  obj
}
