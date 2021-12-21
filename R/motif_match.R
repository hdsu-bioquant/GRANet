
#ArchR:::.summarizeJASPARMotifs()

get_motif_collection <- function(motifSet, genome, collection = "CORE", version = 2){
  species <- list(hg19 = "Homo sapiens",
                  hg38 = "Homo sapiens",
                  mm9  = "Mus musculus",
                  mm10 = "Mus musculus")[[genome]]

  #----------------------------------------------------------------------------#
  #            Get colllection of PWMs - Adapted from ArchR                    #
  #----------------------------------------------------------------------------#
  if(tolower(motifSet)=="jaspar2020"){

    ArchR:::.requirePackage("JASPAR2020",installInfo='BiocManager::install("JASPAR2020")')
    #args <- list(species = species, collection = collection, ...)
    args <- list(species = species, collection = collection)
    motifs <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, args)
    obj <- ArchR:::.summarizeJASPARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary

  }else if(tolower(motifSet)=="jaspar2018"){

    ArchR:::.requirePackage("JASPAR2018",installInfo='BiocManager::install("JASPAR2018")')
    args <- list(species = species, collection = collection)
    motifs <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, args)
    obj <- ArchR:::.summarizeJASPARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary

  }else if(tolower(motifSet)=="jaspar2016"){

    ArchR:::.requirePackage("JASPAR2016",installInfo='BiocManager::install("JASPAR2018")')
    args <- list(species = species, collection = collection)
    motifs <- TFBSTools::getMatrixSet(JASPAR2016::JASPAR2016, args)
    obj <- ArchR:::.summarizeJASPARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary

  }else if(tolower(motifSet)=="cisbp"){

    ArchR:::.requirePackage("chromVARmotifs",installInfo='devtools::install_github("GreenleafLab/chromVARmotifs")')
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
      obj <- ArchR:::.summarizeChromVARMotifs(motifs)
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
      obj <- ArchR:::.summarizeChromVARMotifs(motifs)
      motifs <- obj$motifs
      motifSummary <- obj$motifSummary
    }else{
      stop("Species not recognized homo sapiens, mus musculus supported by CisBP!")
    }

  }else if(tolower(motifSet)=="encode"){

    ArchR:::.requirePackage("chromVARmotifs",installInfo='devtools::install_github("GreenleafLab/chromVARmotifs")')
    data("encode_pwms")
    motifs <- encode_pwms
    obj <- ArchR:::.summarizeChromVARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary

  }else if(tolower(motifSet)=="homer"){

    ArchR:::.requirePackage("chromVARmotifs",installInfo='devtools::install_github("GreenleafLab/chromVARmotifs")')
    data("homer_pwms")
    motifs <- homer_pwms
    obj <- ArchR:::.summarizeChromVARMotifs(motifs)
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
    list(motifs = motifs,
         motifSummary=motifSummary)
  )
}

# y <- get_motif_collection(motifSet = "jaspar2020", genome = "hg19")
# y <- get_motif_collection(motifSet = "encode", genome = "hg19")
# y <- get_motif_collection(motifSet = "homer", genome = "hg19")
# y <- get_motif_collection(motifSet = "cisbp", genome = "hg19")
# y <- get_motif_collection(motifSet = "custom", genome = "hg19")
# y$motifSummary
# names(y)
# table(sub("_.*", "", names(y$motifs)) %in% hs_hgnc_curated_tfs)
# table(y$motifSummary$name  %in% hs_hgnc_curated_tfs)
# table(sub("_.*", "", names(y$motifs)) %in% rownames(granetobj@SeuratObject))
# table(y$motifSummary$name  %in% rownames(granetobj@SeuratObject))


