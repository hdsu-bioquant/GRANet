



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
