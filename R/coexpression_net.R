
#' Co-expression modules with GRNBoost2
#'
#' Computes co-expression modules using GRNBoost2. This function uses the python
#' packages Arboreto and pySCENIC, please install them before (the vignette
#' 01_create_condaenv can be follow to install all the python requirements).
#' The input is a GRANet object initialized from a Seurat object containing the
#' scRNA-seq data as counts, and a list of transcription factors as gene
#' symbols. The gene names in the original Seurat object should be symbols as
#' well.
#'
#' @param GRANetObject Object created using the function CreateGRAnetObject
#' @param threads Number of threads to use.
#'
#' @return GRANetObject with computed co-expression modules.
#' @export
#'
#' @examples
#' \dontrun{
#' granetobj <- coexpression_modules_grnboost2(GRANetObject = granetobj,
#' threads = 8)
#' }
coexpression_modules_grnboost2 <- function(
  GRANetObject,
  threads
){
  GRANetObject <- compute_coexpression_modules(
    GRANetObject = GRANetObject,
    TFs = unique(GRANetObject@TFmotif_location$motifSummary$name),
    threads = threads
  )
  GRANetObject <- add_correlation_to_coexpression_modules(
    GRANetObject,
    mask_dropouts=FALSE
  )

  return(GRANetObject)
}


#' Co-expression modules with GRNBoost2
#'
#' Computes co-expression modules using GRNBoost2. This function uses the python
#' packages Arboreto and pySCENIC, please install them before (the vignette
#' 01_create_condaenv can be follow to install all the python requirements).
#' The input is a GRANet object initialized from a Seurat object containing the
#' scRNA-seq data as counts, and a list of transcription factors as gene
#' symbols. The gene names in the original Seurat object should be symbols as
#' well.
#'
#' @param GRANetObject Object created using the function CreateGRAnetObject
#' @param TFs Character vector with the list of transcription factors included
#' in the Seurat object used to create the GRANetObject.
#' @param threads Number of threads to use.
#'
#' @return GRANetObject with computed co-expression modules.
#' @export
#'
#' @examples
#' \dontrun{
#' granetobj <- compute_coexpression_modules(GRANetObject = granetobj,
#' TFs = readLines("data/mm_mgi_tfs.txt"),
#' threads = 8)
#' }
compute_coexpression_modules <- function(GRANetObject, TFs, threads=1){
  # if(!file.exists(GRANetObject@ProjectMetadata$pathLoom)){
  #   stop("Loom file with counts not found. Did you move the output directory from
  #        its original location? You can reset it with the function
  #        set_GRANet_output_directory")
  # }
  # outputDirectory <- GRANetObject@ProjectMetadata$outputDirectory
  # pathAdjacencies <- file.path(outputDirectory, "Coexprs_modules/expr_mat_adjacencies.tsv")

  # pathLoom <- GRANetObject@ProjectMetadata$pathLoom
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
  cmods <- coexpression_modules(
    method       = "grnboost2",
    ex_matrix_r  = reticulate::r_to_py(Seurat::GetAssayData(
      object = GRANetObject@SeuratObject,
      slot = "counts")),
    gene_names_r = rownames(GRANetObject@SeuratObject),
    tf_names     = TFs,
    num_workers  = as.integer(threads),
    sparse       = TRUE
  )
  GRANetObject@Coexprs_modules <- cmods


  return(GRANetObject)
}




#' Add correlation to co-expression modules
#'
#' Computes the Spearman correlation between the expression of target genes and
#' transcription factors for all co-expression modules.
#' This step is needed to distinguish between positively and negatively
#' regulated co-expression modules.
#'
#' @param GRANetObject GRANet object with computed co-expression modules.
#' @param mask_dropouts Option to include or mask out cells with 0 counts for a
#' given gene when computing the correlations.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' granetobj <- add_correlation_to_coexpression_modules(granetobj, mask_dropouts=FALSE)
#' }
add_correlation_to_coexpression_modules <- function(GRANetObject, mask_dropouts=FALSE){
  if(ncol(GRANetObject@Coexprs_modules) == 0){
    stop("No co-expression modules found. Compute them with the function compute_coexpression_modules")
  }

  if(all(c("TF", "target", "importance", "regulation", "rho") %in% colnames(GRANetObject@Coexprs_modules))){
    stop("Correlation already added to co-expression modules")
  }


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
    ex_matrix = reticulate::r_to_py(Seurat::GetAssayData(
      object = GRANetObject@SeuratObject,
      slot = "counts")),
    gene_names = rownames(GRANetObject@SeuratObject),
    cell_names = colnames(GRANetObject@SeuratObject),
    adj                  = GRANetObject@Coexprs_modules,
    mask_dropouts        = mask_dropouts)
  GRANetObject@Coexprs_modules <- cmods

  return(GRANetObject)
}

# granetobj <- compute_coexpression_modules(GRANetObject = granetobj, TFs = head(readLines("data/mm_mgi_tfs.txt")), threads = 8)
#granetobj <- add_correlation_to_coexpression_modules(granetobj, mask_dropouts=FALSE)
#saveRDS(granetobj, "GRANetProject/granetobj.RDS")

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
