#------------------------------------------------------------------------------#
#            Helper functions  - Taken from ArchR                              #
#------------------------------------------------------------------------------#
.requirePackage <- function (x = NULL, load = TRUE, installInfo = NULL, source = NULL)
{
  if (x %in% rownames(installed.packages())) {
    if (load) {
      suppressPackageStartupMessages(require(x, character.only = TRUE))
    }
    else {
      return(0)
    }
  }
  else {
    if (!is.null(source) & is.null(installInfo)) {
      if (tolower(source) == "cran") {
        installInfo <- paste0("install.packages(\"",
                              x, "\")")
      }
      else if (tolower(source) == "bioc") {
        installInfo <- paste0("BiocManager::install(\"",
                              x, "\")")
      }
      else {
        stop("Unrecognized package source, available are cran/bioc!")
      }
    }
    if (!is.null(installInfo)) {
      stop(paste0("Required package : ", x, " is not installed/found!\n  Package Can Be Installed : ",
                  installInfo))
    }
    else {
      stop(paste0("Required package : ", x, " is not installed/found!"))
    }
  }
}

.summarizeJASPARMotifs <- function (motifs = NULL)
{
  motifNames <- lapply(seq_along(motifs), function(x) {
    namex <- make.names(motifs[[x]]@name)
    if (grepl("LINE", namex)) {
      splitNamex <- stringr::str_split(motifs[[x]]@ID,
                                       pattern = "\\_", simplify = TRUE)
      namex <- splitNamex[1, grep("LINE", splitNamex[1,
      ]) + 1]
    }
    if (substr(namex, nchar(namex), nchar(namex)) == ".") {
      namex <- substr(namex, 1, nchar(namex) - 1)
    }
    namex <- paste0(namex, "_", x)
    namex
  }) %>% unlist(.)
  motifNames2 <- lapply(seq_along(motifs), function(x) {
    namex <- make.names(motifs[[x]]@name)
    if (grepl("LINE", namex)) {
      splitNamex <- stringr::str_split(motifs[[x]]@ID,
                                       pattern = "\\_", simplify = TRUE)
      namex <- splitNamex[1, grep("LINE", splitNamex[1,
      ]) + 1]
    }
    if (substr(namex, nchar(namex), nchar(namex)) == ".") {
      namex <- substr(namex, 1, nchar(namex) - 1)
    }
    namex
  }) %>% unlist(.)
  motifDF <- lapply(seq_along(motifs), function(x) {
    df <- data.frame(row.names = motifNames[x], name = motifNames2[[x]],
                     ID = motifs[[x]]@ID, strand = motifs[[x]]@strand,
                     stringsAsFactors = FALSE)
  }) %>% Reduce("rbind", .) %>% data.frame
  names(motifs) <- motifNames
  out <- list(motifs = motifs, motifSummary = motifDF)
  return(out)
}

.summarizeChromVARMotifs <- function (motifs = NULL)
{
  motifNames <- lapply(seq_along(motifs), function(x) {
    namex <- make.names(motifs[[x]]@name)
    if (grepl("LINE", namex)) {
      splitNamex <- stringr::str_split(motifs[[x]]@ID,
                                       pattern = "\\_", simplify = TRUE)
      namex <- splitNamex[1, grep("LINE", splitNamex[1,
      ]) + 1]
    }
    if (substr(namex, nchar(namex), nchar(namex)) == ".") {
      namex <- substr(namex, 1, nchar(namex) - 1)
    }
    namex <- paste0(namex, "_", x)
    namex
  }) %>% unlist(.)
  motifNames2 <- lapply(seq_along(motifs), function(x) {
    namex <- make.names(motifs[[x]]@name)
    if (grepl("LINE", namex)) {
      splitNamex <- stringr::str_split(motifs[[x]]@ID,
                                       pattern = "\\_", simplify = TRUE)
      namex <- splitNamex[1, grep("LINE", splitNamex[1,
      ]) + 1]
    }
    if (substr(namex, nchar(namex), nchar(namex)) == ".") {
      namex <- substr(namex, 1, nchar(namex) - 1)
    }
    namex
  }) %>% unlist(.)
  motifDF <- lapply(seq_along(motifs), function(x) {
    df <- data.frame(row.names = motifNames[x], name = motifNames2[[x]],
                     ID = motifs[[x]]@ID, strand = motifs[[x]]@strand,
                     stringsAsFactors = FALSE)
  }) %>% Reduce("rbind", .) %>% data.frame
  names(motifs) <- motifNames
  out <- list(motifs = motifs, motifSummary = motifDF)
  return(out)
}


check_BSGenome <- function(
  GRANetObject
){
  genome <- GRANetObject@ProjectMetadata$Genome
  install <- TRUE

  #Check if BSgenome exists!
  if(tolower(genome)=="hg19"){
    if(!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)){
      if(install){
        message("BSgenome for hg19 not installed! Now installing by the following:\n\tBiocManager::install(\"BSgenome.Hsapiens.UCSC.hg19\")")
        BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
      }else{
        stop("BSgenome for hg19 not installed! Please install by setting install = TRUE or by the following:\n\tBiocManager::install(\"BSgenome.Hsapiens.UCSC.hg19\")")
      }
    }
  }else if(tolower(genome)=="hg19test"){
    if(!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)){
      if(install){
        message("BSgenome for hg19 not installed! Now installing by the following:\n\tBiocManager::install(\"BSgenome.Hsapiens.UCSC.hg19\")")
        BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
      }else{
        stop("BSgenome for hg19 not installed! Please install by setting install = TRUE or by the following:\n\tBiocManager::install(\"BSgenome.Hsapiens.UCSC.hg19\")")
      }
    }
  }else if(tolower(genome)=="hg38"){
    if(!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)){
      if(install){
        message("BSgenome for hg38 not installed! Now installing by the following:\n\tBiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\")")
        BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
      }else{
        stop("BSgenome for hg38 not installed! Please install by setting install = TRUE or by the following:\n\tBiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\")")
      }
    }
  }else if(tolower(genome)=="mm9"){
    if(!requireNamespace("BSgenome.Mmusculus.UCSC.mm9", quietly = TRUE)){
      if(install){
        message("BSgenome for mm9 not installed! Now installing by the following:\n\tBiocManager::install(\"BSgenome.Mmusculus.UCSC.mm9\")")
        BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")
      }else{
        stop("BSgenome for mm9 not installed! Please install by setting install = TRUE or by the following:\n\tBiocManager::install(\"BSgenome.Mmusculus.UCSC.mm9\")")
      }
    }
  }else if(tolower(genome)=="mm10"){
    if(!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)){
      if(install){
        message("BSgenome for mm10 not installed! Now installing by the following:\n\tBiocManager::install(\"BSgenome.Mmusculus.UCSC.mm10\")")
        BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
      }else{
        stop("BSgenome for mm10 not installed! Please install by setting install = TRUE or by the following:\n\tBiocManager::install(\"BSgenome.Mmusculus.UCSC.mm10\")")
      }
    }
  }

}
#check_BSGenome(granetobj)


get_motif_collection <- function(motifSet, genome, collection = "CORE", version = 2){
  species <- list(hg19 = "Homo sapiens",
                  hg38 = "Homo sapiens",
                  mm9  = "Mus musculus",
                  mm10 = "Mus musculus")[[genome]]

  #----------------------------------------------------------------------------#
  #            Get colllection of PWMs - Adapted from ArchR                    #
  #----------------------------------------------------------------------------#
  if(tolower(motifSet)=="jaspar2020"){

    .requirePackage("JASPAR2020",installInfo='BiocManager::install("JASPAR2020")')
    #args <- list(species = species, collection = collection, ...)
    args <- list(species = species, collection = collection)
    motifs <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, args)
    obj <- .summarizeJASPARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary

  }else if(tolower(motifSet)=="jaspar2018"){

    .requirePackage("JASPAR2018",installInfo='BiocManager::install("JASPAR2018")')
    args <- list(species = species, collection = collection)
    motifs <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, args)
    obj <- .summarizeJASPARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary

  }else if(tolower(motifSet)=="jaspar2016"){

    .requirePackage("JASPAR2016",installInfo='BiocManager::install("JASPAR2018")')
    args <- list(species = species, collection = collection)
    motifs <- TFBSTools::getMatrixSet(JASPAR2016::JASPAR2016, args)
    obj <- .summarizeJASPARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary

  }else if(tolower(motifSet)=="cisbp"){

    .requirePackage("chromVARmotifs",installInfo='devtools::install_github("GreenleafLab/chromVARmotifs")')
    if(tolower(species) == "mus musculus"){
      if(version == 1){
        message("Using version 1 motifs!")
        data("mouse_pwms_v1")
        motifs <- mouse_pwms_v1
      }else if(version == 2){
        message("Using version 2 motifs!")
        data("mouse_pwms_v2")
        motifs <- mouse_pwms_v2
      }else{
        stop("Only versions 1 and 2 exist!")
      }
      obj <- .summarizeChromVARMotifs(motifs)
      motifs <- obj$motifs
      motifSummary <- obj$motifSummary
    }else if(tolower(species) == "homo sapiens"){
      if(version == 1){
        message("Using version 1 motifs!")
        data("human_pwms_v1")
        motifs <- human_pwms_v1
      }else if(version == 2){
        message("Using version 2 motifs!")
        data("human_pwms_v2")
        motifs <- human_pwms_v2
      }else{
        stop("Only versions 1 and 2 exist!")
      }
      obj <- .summarizeChromVARMotifs(motifs)
      motifs <- obj$motifs
      motifSummary <- obj$motifSummary
    }else{
      stop("Species not recognized homo sapiens, mus musculus supported by CisBP!")
    }

  }else if(tolower(motifSet)=="encode"){

    .requirePackage("chromVARmotifs",installInfo='devtools::install_github("GreenleafLab/chromVARmotifs")')
    data("encode_pwms")
    motifs <- encode_pwms
    obj <- .summarizeChromVARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary

  }else if(tolower(motifSet)=="homer"){

    .requirePackage("chromVARmotifs",installInfo='devtools::install_github("GreenleafLab/chromVARmotifs")')
    data("homer_pwms")
    motifs <- homer_pwms
    obj <- .summarizeChromVARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary

  # }else if(tolower(motifSet)=="custom"){
  #
  #   obj <- NULL
  #   motifs <- motifPWMs
  #   motifSummary <- NULL

  }else{

    stop("Error MotifSet Not Recognized!")

  }


  return(
    list(motifs       = motifs,
         motifSummary = motifSummary)
  )
}

# y <- get_motif_collection(motifSet = "jaspar2020", genome = "hg19")
# y <- get_motif_collection(motifSet = "encode", genome = "hg19")
# y <- get_motif_collection(motifSet = "homer", genome = "hg19")
# y <- GRANet:::get_motif_collection(motifSet = "cisbp", genome = "hg19")
# y <- get_motif_collection(motifSet = "custom", genome = "hg19")
# y$motifSummary
# names(y)
# table(sub("_.*", "", names(y$motifs)) %in% hs_hgnc_curated_tfs)
# table(y$motifSummary$name  %in% hs_hgnc_curated_tfs)
# table(sub("_.*", "", names(y$motifs)) %in% rownames(granetobj@SeuratObject))
# table(y$motifSummary$name  %in% rownames(granetobj@SeuratObject))
