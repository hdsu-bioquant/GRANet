
library(ArchR)


genome <- "hg19"

setwd("/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/")
addArchRGenome(genome)

#------------------------------------------------------------------------------#
#                                find fragments                                #
#------------------------------------------------------------------------------#

# Find all fragments files
fragments_list <- list.files(path = "/media/ag-cherrmann/projects/10_charite_covid19/data/scatac", pattern = "fragments.tsv.gz", recursive = TRUE, full.names = TRUE)
fragments_list <- grep(pattern = "outs/fragments.tsv.gz$", fragments_list, value = TRUE)
# Exclude negative samples
# negative_samples <- c("PBMC-SP5", "SC2-50049", "SC2-50069")
# negative_samples <- c("SC2-50049", "SC2-50069")
negative_samples <- c("PBMC-SP5")
fragments_list <- grep(pattern = paste(negative_samples, collapse="|"), fragments_list, value = TRUE, invert = TRUE)

fragments_list

# Barcode info
barcode_list <- file.path(dirname(fragments_list), "singlecell.csv")
file.exists(barcode_list)
sampleIDs <- basename(sub("/outs/.*", "", barcode_list))



#------------------------------------------------------------------------------#
#                                Build arrow files                             #
#------------------------------------------------------------------------------#
# addArchRThreads(1)
# x <- lapply(13:15, function(i){
#   print(paste("sample: ", i, sampleIDs[i]))
#   
#   ArrowFiles <- createArrowFiles(
#     inputFiles = fragments_list[i],
#     sampleNames = sampleIDs[i],
#     validBarcodes = getValidBarcodes(csvFiles = barcode_list[i],
#                                      sampleNames = sampleIDs[i]),
#     filterTSS = 2, #Dont set this too high because you can always increase later
#     filterFrags = 1000, 
#     addTileMat = TRUE,
#     addGeneScoreMat = TRUE,
#     force = TRUE
#   )
#   ArrowFiles
#   
# })
# 

ArrowFiles <- createArrowFiles(
  inputFiles = fragments_list,
  sampleNames = sampleIDs,
  validBarcodes = getValidBarcodes(csvFiles = barcode_list,
                                   sampleNames = sampleIDs),
  filterTSS = 2, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE
)
ArrowFiles
#ArrowFiles <- c(ArrowFiles, ArrowFiles2)
#ArrowFiles <- list.files(".", pattern = "*.arrow")




# Inferring scATAC-seq Doublets
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)


setwd("~/10_charite_covid19/")

writeLines(ArrowFiles, "/home/bq_aquintero/projects/charite_covid19_TF_activity/data/scatac/ArrowFiles_list.txt")
