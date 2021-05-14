


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
               ArchR            = "list",
               TFmotif_location = "list",
               Coexprs_modules  = "data.frame",
               cssRegulons      = "list",
               SignFeatures = "data.frame" )

)



CreateGRAnetObject <- function(
  SeuratObject,
  cssCluster,
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
  # seuratRNA <- seRNA
  # seuratRNA$Group <- paste0(seRNA@meta.data[,groupRNA])
  # rm(seRNA)
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
  # build_loom(
  #   file.name = p_loomOut,
  #   dgem      = GetAssayData(object = SeuratObject, slot = "counts"),
  #   title     = "scRNA-seq",
  #   #genome    = "hg19",
  #   default.embedding      = seuratobj@reductions$umap@cell.embeddings,
  #   default.embedding.name = "UMAP of full expression matrix"
  # )
  invisible(capture.output(suppressWarnings(SCopeLoomR::build_loom(
    file.name = p_loomOut,
    dgem      = Seurat::GetAssayData(object = SeuratObject, slot = "counts"),
    title     = "scRNA-seq"
  ))))
  SCopeLoomR::close_loom(p_loomOut)

  #saveRDS(seuratobj, p_seuratobjOut)

  GRANetObject <- new(
    Class = 'GRANet',
    SeuratObject = SeuratObject,
    ProjectMetadata = list(outputDirectory=outputDirectory,
                           pathLoom = p_loomOut)
  )

  return(GRANetObject)
}


#' Title
#'
#' @param GRANetObject
#'
#' @return
#' @export
#'
#' @examples
gene_coex_net <- function(GRANetObject){
  if(!file.exists(GRANetObject@ProjectMetadata$pathLoom)){
    stop("Loom file with counts not found. Did you move the output directory from
         its original location? You can reset it with the function
         set_GRANet_output_directory")
  }
  GRANetObject@ProjectMetadata$pathLoom
}
environment(gene_coex_net) <- asNamespace('GRANet')


gene_coex_net(GRANetObject = granetobj)
