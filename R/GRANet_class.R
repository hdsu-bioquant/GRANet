


#' Title
#'
#' @slot SeuratObject Seurat.
#' @slot cssCluster character.
#' @slot ProjectMetadata list.
#' @slot TFmotif_location list.
#' @slot Coexprs_modules data.frame.
#' @slot cssRegulons list.
#' @slot cssRegulonsAUCell matrix.
#' @slot  .
#'
#' @return
#' @export
#'
#' @examples
GRANet <- setClass(
  Class = "GRANet",
  slots = list(SeuratObject      = "Seurat",
               cssCluster        = "character",
               ProjectMetadata   = "list",
               TFmotif_location  = "list",
               Coexprs_modules   = "data.frame",
               cssRegulons       = "list",
               cssRegulonsAUCell = "matrix")

)



CreateGRAnetObject <- function(
  SeuratObject,
  cssCluster,
  genome,
  threads
){

  if(!cssCluster %in% colnames(SeuratObject@meta.data)){
    stop("cssCluster not in meta.data of Seurat Object")
  }

  #------------------------------------#
  #           Add cluster ID           #
  #------------------------------------#
  SeuratObject$cssCluster <- SeuratObject@meta.data[,cssCluster]
  #------------------------------------#
  #           Create loom              #
  #------------------------------------#

  #countsmat <- SeuratObject@assays$RNA@counts
  #annot <- seuratobj@meta.data
  #dim(countsmat)


  GRANetObject <- new(
    Class = 'GRANet',
    SeuratObject = SeuratObject,
    ProjectMetadata = list(Genome = genome)
  )

  return(GRANetObject)
}
#environment(CreateGRAnetObject) <- asNamespace('GRANet')
#granetobj <- CreateGRAnetObject(SeuratObject = seuratobj, cssCluster = "tissue", Force=TRUE )

