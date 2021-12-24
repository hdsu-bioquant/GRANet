


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
    #out     = "positions",
    out     = "matches",
    p.cutoff = cutOff,
    w = width
  )

  GRANetObject@TFmotif_location$motifPositions <- motifPositions
  return(GRANetObject)
  # return(motifPositions)
}



annotate_motif_positions <- function(
  GRANetObject,
  promoter_size
){

  genome <- GRANetObject@ProjectMetadata$Genome

  #----------------------------------------------#
  # Find genes part of the co-expression modules #
  #----------------------------------------------#
  message("Extracting genomic location genes in co-expression modules...")
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

  #1. get gene promoters
  tssGr <- resize(genesGr, width = 1, fix = "start")
  promoters <- trim(suppressWarnings(tssGr+promoter_size))

  # Add peak and target gene
  # Unlist GrangesList and add TF id
  gr    <- GRANetObject@TFmotif_location$motifPositions
  idx   <- match(names(gr), rownames(GRANetObject@TFmotif_location$motifSummary))
  gr$TF <- GRANetObject@TFmotif_location$motifSummary$name[idx]
  # Overlap with promoters and keep only binding motifs around TSS
  grhits <- findOverlaps(query = gr, subject = promoters)
  gr_promoters <- gr[queryHits(grhits)]
  gr_promoters$Target <- promoters$symbol[subjectHits(grhits)]
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
