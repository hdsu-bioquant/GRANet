#' Find TFs motif positions linked to open chromatin regions and make regulons
#'
#' This function takes as input a GRANetMultiome object,
#' and it extracts the position of motifs associated to transcription factors
#' linked to open chromatin regions.
#' All the positions of transcription factors motifs included in the co-expression
#' modules stored in the GRANetMultiome object are searched using motifmatchr.
#' Then the cis-regulatory interactions inside the selected promoter size, are
#' kept in order to create the regulons.
#'
#'
#' @param GRANetObject GRANet object with computed co-expression modules and
#' Spearman correlation between transcription factors and target genes.
#' @param promoter_size Window size to search for transcription factor motifs
#' located around the TSS of all target genes from the co-expression modules.
#' @param min_regulon_size Minimum number of target genes required to build a
#' regulon, composed of a transcription factor and a list of cis-regulated
#' target genes.
#' @param cutOff motifmatchr.: p-value cutoff for returning motifs
#' @param width motifmatchr: p-value cutoff for returning motifs

#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' granetobj <- match_motifs_and_trim(GRANetObject=granetobj,
#' promoter_size=10000)
#' }
match_motifs_and_trim <- function(
  GRANetObject,
  promoter_size,
  min_regulon_size=20,
  cutOff = 5e-05,
  width = 7
){
  message("Finding position of motifs in open chromatin sites")
  GRANetObject <- find_motif_positions(GRANetObject = GRANetObject,
                                       cutOff       = cutOff,
                                       width        = width)
  message("Annotating motifs positions with target genes and peaksID")
  GRANetObject <- annotate_motif_positions(GRANetObject  = GRANetObject,
                                           promoter_size = promoter_size)

  message("Building regulons")
  GRANetObject <- make_regulons(GRANetObject     = GRANetObject,
                                min_regulon_size = min_regulon_size)

  return(GRANetObject)
}


find_motif_positions <- function(
  GRANetObject,
  cutOff = 5e-05,
  width = 7
){
  genome <- GRANetObject@ProjectMetadata$Genome
  check_BSGenome(GRANetObject)

  #----------------------------------------------#
  #            Load genome sequence              #
  #----------------------------------------------#
  if (genome == "hg19") {
    library(BSgenome.Hsapiens.UCSC.hg19)
    BSgenome <- BSgenome.Hsapiens.UCSC.hg19
  } else if (genome == "hg38") {
    library(BSgenome.Hsapiens.UCSC.hg38)
    BSgenome <- BSgenome.Hsapiens.UCSC.hg38
  } else if (genome == "mm9") {
    library(BSgenome.Mmusculus.UCSC.mm9)
    BSgenome <- BSgenome.Mmusculus.UCSC.mm9
  } else if (genome == "mm10") {
    library(BSgenome.Mmusculus.UCSC.mm10)
    BSgenome <- BSgenome.Mmusculus.UCSC.mm10
  }

  motifPositions <- motifmatchr::matchMotifs(
    pwms    = GRANetObject@TFmotif_location$motifs,
    subject = GRANetObject@TFmotif_location$peaks_GRanges,
    genome  = BSgenome,
    out     = "positions",
    #out     = "matches",
    p.cutoff = cutOff,
    w = width
  )

  GRANetObject@TFmotif_location$motifPositions <- motifPositions
  return(GRANetObject)
  # return(motifPositions)
}
environment(find_motif_positions) <- asNamespace("GRANet")


annotate_motif_positions <- function(
  GRANetObject,
  promoter_size
){

  genome <- GRANetObject@ProjectMetadata$Genome

  #----------------------------------------------#
  # Find genes part of the co-expression modules #
  #----------------------------------------------#
  message("Extracting genomic location for genes in co-expression modules...")
  if (genome == "hg19") {
    library(EnsDb.Hsapiens.v75)
    genesGr <- genes(EnsDb.Hsapiens.v75)
  } else if (genome == "hg38") {
    library(EnsDb.Hsapiens.v86)
    genesGr <- genes(EnsDb.Hsapiens.v86)
  } else if (genome == "mm9") {
    library(EnsDb.Mmusculus.v75)
    genesGr <- genes(EnsDb.Mmusculus.v75)
  } else if (genome == "mm10") {
    library(EnsDb.Mmusculus.v79)
    genesGr <- genes(EnsDb.Mmusculus.v79)
  }
  seqlevelsStyle(genesGr) <- 'UCSC'

  # Get gene promoters
  message("Add TF ID...")
  tssGr <- resize(genesGr, width = 1, fix = "start")
  promoters <- trim(suppressWarnings(tssGr+promoter_size))

  # Add peak and target gene
  # Unlist GrangesList and add TF id
  gr    <- unlist(GRANetObject@TFmotif_location$motifPositions)
  idx   <- match(names(gr), rownames(GRANetObject@TFmotif_location$motifSummary))
  gr$TF <- GRANetObject@TFmotif_location$motifSummary$name[idx]

  # Overlap with promoters and keep only binding motifs around TSS
  message("Add peak and target gene...")
  grhits <- findOverlaps(query = gr, subject = promoters)
  gr_promoters <- gr[queryHits(grhits)]
  gr_promoters$target <- promoters$symbol[subjectHits(grhits)]
  names(gr_promoters) <- NULL

  # Overlap with peaks to find peak id
  #peakIDs <- paste0("P", 1:length(GRANetObject@TFmotif_location$peaks_GRanges))
  #GRANetObject@TFmotif_location$peaks_GRanges$PeakID <- peakIDs
  grhits <- findOverlaps(query   = gr_promoters,
                         subject = GRANetObject@TFmotif_location$peaks_GRanges)
  gr_prom_peaks <- gr_promoters[queryHits(grhits)]
  gr_prom_peaks$PeakID <- GRANetObject@TFmotif_location$peaks_GRanges$PeakID[subjectHits(grhits)]


  GRANetObject@TFmotif_location$annotated_motifs <- gr_prom_peaks
  return(GRANetObject)

}

make_regulons <- function(
  GRANetObject,
  min_regulon_size=20
){

  # TF-target cis-regulatory interactions
  TF_target_motifs <- GRANetObject@TFmotif_location$annotated_motifs %>%
    as.data.frame() %>%
    dplyr::select(c("TF", "target")) %>%
    dplyr::distinct()


  regulons <- GRANetObject@Coexprs_modules %>%
    # keep only positve and negative correlation
    dplyr::filter(!regulation == 0) %>%
    # keep only TF-target cis-regulatory interactions
    inner_join(TF_target_motifs, by = c("TF", "target")) %>%
    # Remove small regulons
    dplyr::group_by(TF, regulation) %>%
    dplyr::filter(n() >= min_regulon_size) %>%
    # Split by regulon
    dplyr::mutate(RegulonID = paste0(TF, " (",
                                     if_else(regulation == 1, '+', "-"),
                                     ")")) %>%
    named_group_split(RegulonID) %>%
    # add TF to regulon
    purrr::map(function(x) c(unique(x$TF), x$target))

  GRANetObject@Regulons <- regulons

  return(GRANetObject)
}

