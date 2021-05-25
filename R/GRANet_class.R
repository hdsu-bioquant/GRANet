


#' Creates a new GRANet object
#'
#' @slot SeuratObject Seurat.
#' @slot cssCluster character.
#' @slot ProjectMetadata list.
#' @slot TFmotif_location list.
#' @slot Coexprs_modules data.frame.
#' @slot cssRegulons list.
#' @slot cssRegulonsAUCell matrix.
#'
#' @importClassesFrom Seurat Seurat
#'
#' @return
#' @export
#'
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



#' Creates a GRAnet object
#'
#' A GRANet object is initialized from a Seurat object containing scRNA-seq
#' data as counts. The gene names in the original count matrix should be gene
#' symbols as they are required to match to the motifs associated to a list of
#' transcription factors. Additionally, one of the columns of the cell metadata
#' from the Seurat object should contain a categorical annotation of the cells
#' (e.g., cluster, cell type, cell state).
#'
#' @param SeuratObject Seurat object with scRNA-seq counts and a categorical
#' cell identity annotation.
#' @param cssCluster Name of the column in the Seurat object that contains the
#' categorical cell identity annotation (e.g., cluster, cell type, cell state).
#' @param genome Genome of the organism under study, supported genomes are:
#' "hg19", "hg38", "mm9", and "mm10".
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' granetobj <- CreateGRAnetObject(SeuratObject = seuratobj,
#' cssCluster = "tissue",
#' genome = "mm9")
#' }
CreateGRAnetObject <- function(
  SeuratObject,
  cssCluster,
  genome
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

