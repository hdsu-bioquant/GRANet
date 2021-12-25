


#' Creates a new GRANet object
#'
#' @slot SeuratObject Seurat.
#' @slot ProjectMetadata list.
#' @slot TFmotif_location list.
#' @slot Coexprs_modules data.frame.
#' @slot cssRegulons list.
#' @slot cssRegulonsAUCell matrix.
#'
#' @importClassesFrom Seurat Seurat
#'
#' @return
#' @export
#'
GRANetMultiome <- setClass(
  Class = "GRANetMultiome",
  slots = list(SeuratObject     = "Seurat",
               ProjectMetadata  = "list",
               TFmotif_location = "list",
               Coexprs_modules  = "data.frame",
               Regulons         = "list",
               RegulonsAUCell   = "list",
               Reg_interactions = "list")
)



#' Creates a GRAnet object
#'
#' A GRANetMultiome object is initialized from a Seurat object containing scRNA-seq
#' data as counts. The gene names in the original count matrix should be gene
#' symbols as they are required to match to the motifs associated to a list of
#' transcription factors. Additionally, one of the columns of the cell metadata
#' from the Seurat object should contain a categorical annotation of the cells
#' (e.g., cluster, cell type, cell state).
#'
#' @param SeuratObject Seurat object with scRNA-seq counts and a categorical
#' cell identity annotation.
#' @param cssCluster Name of the column in the Seurat object that contains the
#' categorical cell identity annotation (e.g., cluster, cell type, cell state).
#' @param genome Genome of the organism under study, supported genomes are:
#' "hg19", "hg38", "mm9", and "mm10".
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' granetobj <- CreateGRAnetObject(SeuratObject = seuratobj,
#' cssCluster = "tissue",
#' genome = "mm9")
#' }
CreateGRANetMultiomeObject <- function(
  SeuratObject,
  peak_counts,
  peaks_GRanges,
  motifSet,
  genome
){
  # message("Extracting genomic location genes in co-expression modules...")
  if (!genome %in% c("hg19", "hg38", "mm9", "mm10")) {
    stop("Non supported genome, please use one of the following:\n",
         "hg19, hg38, mm9, mm10")
  }
  if (is.null(peaks_GRanges$PeakID)) {
    stop("PeakID column missing in peaks_GRanges",
         "This column should match exactly the row names on the peak_counts matrix")
  }
  if (!identical(peaks_GRanges$PeakID, rownames(peak_counts))) {
    stop("PeakID column in peaks_GRanges does not match peaks_GRanges row names ",
         "This column should match exactly the row names on the peak_counts matrix")
  }



  # if(!cssCluster %in% colnames(SeuratObject@meta.data)){
  #   stop("cssCluster not in meta.data of Seurat Object")
  # }

  # get TF motifs and keep only those included in the scRNA-seq data
  motifs <- get_motif_collection(motifSet = motifSet, genome = "hg19")
  message("A total of ", length(motifs$motifs), " motifs found in ", motifSet)

  idx <- motifs$motifSummary$name %in% rownames(SeuratObject)
  motifs$motifSummary <- motifs$motifSummary[idx,]
  idx <- names(motifs$motifs) %in% rownames(motifs$motifSummary)
  motifs$motifs <- motifs$motifs[idx]
  message("Using only ", length(motifs$motifs),
          " motifs for which a corresponding TF was found in the scRNA-seq data")

  # transform peak matrix to binary sparse matrix
  peak_counts <- as(peak_counts, "lgCMatrix")


  GRANetObject <- new(
    Class            = 'GRANetMultiome',
    SeuratObject     = SeuratObject,
    TFmotif_location = list(motifs        = motifs$motifs,
                            motifSummary  = motifs$motifSummary,
                            # Add Granges and peaks counts to GRANet object
                            peak_counts   = peak_counts,
                            peaks_GRanges = peaks_GRanges
                            ),
    ProjectMetadata  = list(Genome = genome)
  )

  return(GRANetObject)
}
#environment(CreateGRAnetObject) <- asNamespace('GRANet')
#granetobj <- CreateGRAnetObject(SeuratObject = seuratobj, cssCluster = "tissue", Force=TRUE )

