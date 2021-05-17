


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
  #SCopeLoomR::close_loom(p_loomOut)
  #SCopeLoomR::close_loom(local_loom)
  SCopeLoomR::finalize(local_loom)
  print(local_loom)
  #SCopeLoomR::flush(local_loom)
  #SCopeLoomR::close_loom(local_loom)

  #saveRDS(seuratobj, p_seuratobjOut)

  GRANetObject <- new(
    Class = 'GRANet',
    SeuratObject = SeuratObject,
    ProjectMetadata = list(outputDirectory=outputDirectory,
                           pathLoom = p_loomOut)
  )

  return(GRANetObject)
}
environment(CreateGRAnetObject) <- asNamespace('GRANet')
granetobj <- CreateGRAnetObject(SeuratObject = seuratobj, cssCluster = "tissue", Force=TRUE )
#' Title
#'
#' @param GRANetObject
#'
#' @return
#' @export
#'
#' @examples
gene_coex_net <- function(GRANetObject, TFs, threads=1){
  if(!file.exists(GRANetObject@ProjectMetadata$pathLoom)){
    stop("Loom file with counts not found. Did you move the output directory from
         its original location? You can reset it with the function
         set_GRANet_output_directory")
  }
  outputDirectory <- GRANetObject@ProjectMetadata$outputDirectory
  pathAdjacencies <- file.path(outputDirectory, "Coexprs_modules/expr_mat_adjacencies.tsv")

  pathLoom <- GRANetObject@ProjectMetadata$pathLoom
  gnrboos2 <- arboreto_with_multiprocessing.py

  pathTFs

  cmd <- sprintf("-o %s --num_workers %s %s %s",
                 pathAdjacencies, threads, pathLoom, pathTFs)

  "arboreto_with_multiprocessing.py -o /home/bq_aquintero/projects/charite_covid19_TF_activity/atac17_rna57/results/scrna/SCENIC/expr_mat.adjacencies.tsv --num_workers 54 /home/bq_aquintero/projects/charite_covid19_TF_activity/atac17_rna57/results/scrna/SCENIC/rnaseq_counts.loom /home/bq_aquintero/projects/GRANet/src/GRANet/scGRANet/aux/scenic/allTFs_hg38.txt"


  run <- system2(pathToMacs2, cmd, wait=TRUE, stdout=NULL, stderr=NULL)

}
#
# "GRANetProject/Coexprs_modules/mm_mgi_tfs.txt"
#
# reticulate::source_python("inst/arboreto_with_multiprocessing.py -o GRANetProject/Coexprs_modules/expr_mat.adjacencies.tsv --num_workers 1 GRANetProject/Coexprs_modules/counts.loom GRANetProject/Coexprs_modules/mm_mgi_tfs.txt", envir=NULL)
#
# system('source ~/.bashrc && conda activate pyscenic && conda info -e')
#
# system("bash -c 'conda activate pyscenic; conda info -e'")
# system("bash -c 'conda info -e'")
#
# system2(pathToMacs2, cmd, wait=TRUE, stdout=NULL, stderr=NULL)
# system2("arboreto_with_multiprocessing.py")
#
# environment(gene_coex_net) <- asNamespace('GRANet')
# reticulate::miniconda_path()
# compute_coexpression_modules
#
# gene_coex_net(GRANetObject = granetobj)
#
# reticulate::py_config()
