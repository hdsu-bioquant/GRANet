


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
GRANet <- setClass(
  Class = "GRANet",
  slots = list(GeneExpression   = "list",
               SeuratObject     = "Seurat",
               cssCluster       = "character",
               ProjectMetadata  = "list",
               TFmotif_location = "list",
               Coexprs_modules  = "data.frame",
               cssRegulons      = "list")

)



CreateGRAnetObject <- function(
  SeuratObject,
  cssCluster,
  genome,
  outputDirectory = "GRANetProject",
  Force=FALSE,
  threads
){
  if(dir.exists(outputDirectory) & Force==FALSE){
    stop("GRANet output directory already exists, use Force=TRUE to overwrite")
  }

  # Create output directory
  dir.create(outputDirectory, showWarnings=FALSE, recursive=TRUE)


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

  p_loomOut <- file.path(outputDirectory, "Coexprs_modules/counts.loom")
  if(file.exists(p_loomOut)){
    unlink(dirname(p_loomOut), recursive = TRUE)
  }

  dir.create(dirname(p_loomOut), showWarnings=FALSE, recursive=TRUE)

  ### Create the minimal loom file
  # local_loom <- invisible(capture.output(suppressWarnings(SCopeLoomR::build_loom(
  #   file.name = p_loomOut,
  #   dgem      = Seurat::GetAssayData(object = SeuratObject, slot = "counts"),
  #   title     = "scRNA-seq"
  # ))))
  local_loom <- invisible(suppressWarnings(SCopeLoomR::build_loom(
    file.name = p_loomOut,
    dgem      = Seurat::GetAssayData(object = SeuratObject, slot = "counts"),
    title     = "scRNA-seq"
  )))
  SCopeLoomR::finalize(local_loom)

  GRANetObject <- new(
    Class = 'GRANet',
    SeuratObject = SeuratObject,
    ProjectMetadata = list(outputDirectory = outputDirectory,
                           pathLoom        = p_loomOut,
                           Genome          = genome)
  )

  return(GRANetObject)
}
#environment(CreateGRAnetObject) <- asNamespace('GRANet')
#granetobj <- CreateGRAnetObject(SeuratObject = seuratobj, cssCluster = "tissue", Force=TRUE )

