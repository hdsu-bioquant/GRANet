


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
#environment(CreateGRAnetObject) <- asNamespace('GRANet')
#granetobj <- CreateGRAnetObject(SeuratObject = seuratobj, cssCluster = "tissue", Force=TRUE )

#' Title
#'
#' @param GRANetObject
#'
#' @return
#' @export
#'
#' @examples
compute_coexpression_modules <- function(GRANetObject, TFs, threads=1){
  if(!file.exists(GRANetObject@ProjectMetadata$pathLoom)){
    stop("Loom file with counts not found. Did you move the output directory from
         its original location? You can reset it with the function
         set_GRANet_output_directory")
  }
  outputDirectory <- GRANetObject@ProjectMetadata$outputDirectory
  pathAdjacencies <- file.path(outputDirectory, "Coexprs_modules/expr_mat_adjacencies.tsv")

  pathLoom <- GRANetObject@ProjectMetadata$pathLoom
  #gnrboos2 <- arboreto_with_multiprocessing.py

  #pathTFs


  #------------------------------------#
  #   Import GRNBoost2 function        #
  #------------------------------------#
  message("Loading GRNBoost2...")
  path <- system.file(package = "GRANet")
  reticulatedBoost <- reticulate::import_from_path("reticulatedBoost", path = path)
  coexpression_modules <- reticulatedBoost$arboreto_functionalized$coexpression_modules

  #------------------------------------#
  #   Compute coexpression modules     #
  #------------------------------------#
  message("Computing Co-expression modules...")
  cmods <- coexpression_modules(method               = "grnboost2",
                                expression_mtx_fname = pathLoom,
                                tf_names             = TFs,
                                num_workers          = as.integer(threads),
                                sparse               = TRUE)
  GRANetObject@Coexprs_modules <- cmods


  return(GRANetObject)
}




#' Title
#'
#' @param GRANetObject
#' @param threads
#'
#' @return
#' @export
#'
#' @examples
add_correlation_to_coexpression_modules <- function(GRANetObject, mask_dropouts=FALSE){
  if(!file.exists(GRANetObject@ProjectMetadata$pathLoom)){
    stop("Loom file with counts not found. Did you move the output directory from
         its original location? You can reset it with the function
         set_GRANet_output_directory")
  }
  if(ncol(GRANetObject@Coexprs_modules) == 0){
    stop("No co-expression modules found. Compute them with the function compute_coexpression_modules")
  }

  if(all(c("TF", "target", "importance", "regulation", "rho") %in% colnames(GRANetObject@Coexprs_modules))){
    stop("Correlation already added to co-expression modules")
  }


  pathLoom <- GRANetObject@ProjectMetadata$pathLoom

  #------------------------------------#
  #   Import correlation function      #
  #------------------------------------#
  path <- system.file(package = "GRANet")
  reticulatedBoost <- reticulate::import_from_path("reticulatedBoost", path = path)
  correlations_to_modules <- reticulatedBoost$arboreto_functionalized$correlations_to_modules

  #------------------------------------#
  #   Compute coexpression modules     #
  #------------------------------------#
  message("Adding correlation to co-expression modules...")
  cmods <- correlations_to_modules(
    expression_mtx_fname = pathLoom,
    adj                  = GRANetObject@Coexprs_modules,
    mask_dropouts        = mask_dropouts)
  GRANetObject@Coexprs_modules <- cmods

  return(GRANetObject)
}

# granetobj <- compute_coexpression_modules(GRANetObject = granetobj, TFs = head(readLines("data/mm_mgi_tfs.txt")), threads = 8)
# add_correlation_to_coexpression_modules(granetobj, mask_dropouts=FALSE)

#environment(compute_coexpression_modules) <- asNamespace('GRANet')
#compute_coexpression_modules(GRANetObject = granetobj, TFs = readLines("data/mm_mgi_tfs.txt"))
#compute_coexpression_modules(GRANetObject = granetobj, TFs = head(readLines("data/mm_mgi_tfs.txt")), threads = 8)
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
