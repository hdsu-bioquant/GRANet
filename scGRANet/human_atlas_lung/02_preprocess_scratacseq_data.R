#options(echo=TRUE) # if you want see commands in output file
#args <- commandArgs(TRUE)

# p_scenic_adjacencies <- as.character(args[1])

#------------------------------------------------------------------------------#
params <- list(
  genome      = "hg38",
  rawFilesdir = "/media/ag-cherrmann/projects/10_charite_covid19/subprojects/human_atlas_TF_activity/data/atac/fragments/",
  rawFiles    = c("lung_SM-A62E9_rep1_fragments.bed.gz", "lung_SM-A8WNH_rep1_fragments.bed.gz", 
                  "lung_SM-ACCPU_rep1_fragments.bed.gz", "lung_SM-JF1NZ_rep1_fragments.bed.gz"),
  # rawFiles    = c("lung_SM-A62E9_rep1_fragments.txt", "lung_SM-A8WNH_rep1_fragments.txt", 
  #                 "lung_SM-ACCPU_rep1_fragments.txt", "lung_SM-JF1NZ_rep1_fragments.txt"),
  archerOutdir = "/media/ag-cherrmann/projects/10_charite_covid19/subprojects/human_atlas_TF_activity/data/atac/fragments_processed/"
               )

rawFiles <- c("/media/ag-cherrmann/projects/10_charite_covid19/subprojects/human_atlas_TF_activity/data/atac/fragments_processed/lung_SM-A62E9_rep1_fragments.txt",
              "/media/ag-cherrmann/projects/10_charite_covid19/subprojects/human_atlas_TF_activity/data/atac/fragments_processed/lung_SM-A8WNH_rep1_fragments.txt",
              "/media/ag-cherrmann/projects/10_charite_covid19/subprojects/human_atlas_TF_activity/data/atac/fragments_processed/lung_SM-ACCPU_rep1_fragments.txt",
              "/media/ag-cherrmann/projects/10_charite_covid19/subprojects/human_atlas_TF_activity/data/atac/fragments_processed/lung_SM-JF1NZ_rep1_fragments.txt")


archerOutdir <- unique(dirname(rawFiles))
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                               Index and sort fragments                       #
#------------------------------------------------------------------------------#

# Original files are not sorted, and bgzipped
# first read them into a granges, sort, and then index

makeFragTabix<-function(filepath,skip=0,rm.file=TRUE){
  message("compressing the file with bgzip...")
  zipped <- Rsamtools::bgzip(filepath,overwrite=TRUE)
  
  if(rm.file){file.remove(filepath)}
  
  message("making tabix index...")
  Rsamtools::indexTabix(zipped,
                        seq=1, start=2, end=3,
                        skip=skip, comment="#", zeroBased=FALSE)
  
}

rawFiles <- lapply(params$rawFiles, function(rawFile){
  message("Start file: ", rawFile)
  
  from <-  file.path(params$rawFilesdir, rawFile)
  baseto <- paste0(sub("\\..*", "", rawFile), ".txt.gz")
  to <- file.path(params$archerOutdir, baseto)
  
  # Read data
  x <- read.table(from)
  x <- x[,1:5]
  colnames(x) <- c("chr", "start", "end", "barcode", "count")
  # Make granges and sort
  gr <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  gr <-sort(gr, ignore.strand=TRUE)
  # Save as table
  message("save sorted granges as table: ", rawFile)
  y <- as.data.frame(gr)
  y <- y[,c(1,2,3,6,7)]
  write_delim(y, file=to, delim = "\t", col_names = FALSE)
  
  # bgzip and index
  message("bgzip and index: ", rawFile)
  makeFragTabix(to)
  
  # final file
  baseto <- paste0(sub("\\..*", "", rawFile), ".txt.bgz")
  file.path(params$archerOutdir, baseto)
  
})







#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#

library(ArchR)
# library(patchwork)
# library(ComplexHeatmap)
# library(viridis)

addArchRGenome(genome)




myrealwd <- getwd()
setwd(archerOutdir)


# Name sample name from file
sampleNames <- basename(rawFiles)

# sample_ids_atac <- c("SC2-10025AT-NS02FU", "SC2-20017AT-NS02FU", "SC2-30011AT-NS02FU")
# names(sample_ids_atac) <- c("SC2-10025-NS02FU",   "SC2-30011-NS02FU", "SC2-20017-NS02FU")
#getValidBarcodes

x <- read.table("/media/ag-cherrmann/projects/10_charite_covid19/subprojects/human_atlas_TF_activity/data/atac/fragments_processed/lung_SM-A62E9_rep1_fragments.bed.gz")
head(x)
x <- x[,1:5]
colnames(x) <- c("chr", "start", "end", "barcode", "count")
gr <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
gr <-sort(gr, ignore.strand=TRUE)
gr
y <- as.data.frame(gr)
head(y)
y <- y[,c(1,2,3,6,7)]


makeFragTabix<-function(filepath,skip=0,rm.file=TRUE){
  message("compressing the file with bgzip...")
  zipped <- Rsamtools::bgzip(filepath,overwrite=TRUE)
  
  #if(rm.file){file.remove(filepath)}
  
  message("making tabix index...")
  Rsamtools::indexTabix(zipped,
                        seq=1, start=2, end=3,
                        skip=skip, comment="#", zeroBased=FALSE)
  
}

makeFragTabix("/media/ag-cherrmann/projects/10_charite_covid19/subprojects/human_atlas_TF_activity/data/atac/fragments_processed/lung_SM-A62E9_rep1_fragments.txt.gz")

write_delim(y, file=out_tabix, delim = "\t", col_names = FALSE)
indexTabix(out_tabix, seq=1, start=2, end=3)

Rsamtools::bgzip()


out_tabix <- "/media/ag-cherrmann/projects/10_charite_covid19/subprojects/human_atlas_TF_activity/data/atac/fragments_processed/lung_SM-A62E9_rep1_fragments.txt.bgz"
#exportToTabix(gr, con=out_tabix,  quote=FALSE)

#xt <- TabixFile(file = out_tabix)
#indexTabix(out_tabix, format = "bed")

ArrowFiles <- createArrowFiles(
  inputFiles = out_tabix,
  sampleNames = basename(out_tabix),
  minTSS   = 2, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE,
  threads = 1
)
ArrowFiles



ArrowFiles <- createArrowFiles(
  inputFiles = rawFiles,
  sampleNames = sampleNames,
  minTSS   = 2, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE,
  threads = 1
)
ArrowFiles
#ArrowFiles <- list.files(".", pattern = "*.arrow")



# Inferring scATAC-seq Doublets
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)


#setwd("~/10_charite_covid19/")
setwd(myrealwd)
