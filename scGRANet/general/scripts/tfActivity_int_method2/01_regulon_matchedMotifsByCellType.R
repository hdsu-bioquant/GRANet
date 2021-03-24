options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

p_scenic_adjacencies <- as.character(args[1])
p_archrproj          <- dirname(as.character(args[2]))
p_regulons           <- as.character(args[3])
genome               <- as.character(args[4])
importance_threshold <- as.numeric(args[5])/100
promoter_size        <- as.numeric(args[6])
min_regulon_size     <- as.numeric(args[7])
Ncores               <- as.numeric(args[8])
cellWithPeak         <- as.numeric(args[9])
#------------------------------------------------------------------------------#
# p_scenic_adjacencies <- "results/scrna/SCENIC/expr_mat.adjacencies.cor.tsv"
# p_archrproj          <- "results/scatac/archr/Save-projcovid6/"
# #p_archrproj <- "~/projects/mouse_atlas_TF_activity/adult_lung/results/scatac/archr/ArchR02_MotifMatch/"
# #p_scenic_corrmodules <- "results/integrated/TF_activity_method2/tfModules_asDF.RDS"
# p_regulons           <- "results/integrated/TF_activity_method2/tfRegulons_asDF_byCell.RDS"
# genome               <- "hg19"
# importance_threshold <- 50/100
# promoter_size        <- 5000
# min_regulon_size     <- 10
# Ncores               <- 20
# cellWithPeak         <- 0.1 # Minimum fraction of cells for a given cell type, should have this peak to conserve it


# p_scenic_adjacencies <- '~/projects/mouse_atlas_TF_activity/adult_9tissues/results/scrna/SCENIC/expr_mat.adjacencies.cor.tsv'
# p_archrproj          <- '~/projects/mouse_atlas_TF_activity/adult_9tissues/results/scatac/archr/ArchR02_MotifMatch/'
# p_regulons           <- '~/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_5000/tfRegulons_asDF.RDS'
# genome               <- "mm9"
# importance_threshold <- 50/100
# promoter_size        <- 5000
# min_regulon_size     <- 20
# #Ncores               <- as.numeric(args[8])
# cellWithPeak         <- 0.1

#p_archrproj <- "~/projects/mouse_atlas_TF_activity/adult_9tissues/results/scatac/archr/ArchR02_MotifMatch/"
#------------------------------------------------------------------------------#
library(tidyverse)
library(ArchR)

#addMotifAnnotations()
# change cutoff to 0.05

#------------------------------------------------------------------------------#
#               Read adjacencies and make co-expression modules                #
#------------------------------------------------------------------------------#
scenic_adjacencies <- read_tsv(p_scenic_adjacencies)

scenic_corrmodules <- scenic_adjacencies %>% 
  # remove targets with low correlation
  dplyr::filter(!regulation == 0) %>% 
  # Split For each TF and in positive or negative regulation
  group_by(TF, regulation) %>% 
  # Keep only targets with top x% importance
  dplyr::filter(importance >= quantile(importance, probs = importance_threshold))

# Summary
scenic_corrmodules %>% 
  summarise(n = n()) %>% 
  summary()

# Save co-expression modules from arboreto
#saveRDS(scenic_corrmodules, p_scenic_corrmodules)

#------------------------------------------------------------------------------#
#                   Find genes part of the co-expression modules               #
#------------------------------------------------------------------------------#
# Load gr with genes
# addArchRGenome(genome)
# genesGr <- ArchR::getGenes()
# genesGr

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
#genesGr

table(toupper(genesGr$symbol) %in% toupper(unique(scenic_corrmodules$target)))
table(genesGr$symbol %in% unique(scenic_corrmodules$target))

# Find genes in co-expression modules
genesGr <- genesGr[genesGr$symbol %in% unique(scenic_corrmodules$target)]
#genesGr <- genesGr[toupper(genesGr$symbol) %in% toupper(unique(scenic_corrmodules$target))]

# get TSS
tssGr <- resize(genesGr, width = 1, fix = "start")
seqlevelsStyle(tssGr) <- 'UCSC'
tssGr

#------------------------------------------------------------------------------#
#           Filter TFs that have no motif in the ATACseq data                  #
#------------------------------------------------------------------------------#
archrproj <- loadArchRProject(path = p_archrproj)
motifPositions <- getPositions(archrproj)
motifPositions


#getPositions(archrproj3)
#names(motifPositions)

scenicTFs <- tibble(motif_archr = names(motifPositions),
                    TF_archr    = sub("_.*", "", names(motifPositions))) %>% 
  dplyr::mutate(TF = scenic_corrmodules$TF[match(TF_archr, scenic_corrmodules$TF)]) %>% 
  dplyr::filter(!is.na(TF)) %>% 
  dplyr::distinct() 

scenicTFs
# Remove modules where the TF was not found in the ATACseq data
scenic_corrmodules <- scenic_corrmodules %>% 
  dplyr::filter(TF %in% scenicTFs$TF)

# Remove motifs for TF not found 
length(motifPositions)
motifPositions <- motifPositions[names(motifPositions) %in% scenicTFs$motif_archr]
length(motifPositions)
#------------------------------------------------------------------------------#
#                 Find genes with a motif in its promoter.                     #
#------------------------------------------------------------------------------#
# Find peaks between 5k of tss
# To identify genes that are cis regulated by the TF
# And later remove genes that are downstream target of the TF

# Split peaks into peaks found by Cell type
# This way cell type-specific regulos will be buils


# Load peaks for each cell type
peaksGr <- ArchR::getPeakSet(archrproj)
unique(names(peaksGr))

# Get peaks matrix
ArchR::getAvailableMatrices(archrproj)
peakMatrix <- ArchR::getMatrixFromProject(archrproj, useMatrix = "PeakMatrix")
peakMatrix <- assay(peakMatrix)
dim(peakMatrix)
length(peaksGr)

# Get peaks by cell type
#cellWithPeak <- 0.01
table(archrproj$predictedGroup)
cellTypes <- sort(unique(archrproj$predictedGroup))
names(cellTypes) <- cellTypes
peaksGr_CellTypes <- lapply(cellTypes, function(cellType){
  # Keep peak if at least one insertion in x% of cells from group
  peakMatrixCellType <- peakMatrix[,archrproj$predictedGroup == cellType]
  percentCells         <- ncol(peakMatrixCellType)*cellWithPeak
  nCellWithPeak      <- rowSums(peakMatrixCellType > 0)
  peaksGr[nCellWithPeak > percentCells]
})

# Create list of motif position by cell type
motifPositions_Cell <- lapply(peaksGr_CellTypes, function(peaksGr_Cell){
  endoapply(motifPositions, subsetByOverlaps, peaksGr_Cell)
})
#motifPositions_Cell
sapply(motifPositions_Cell, function(x){
  summary(sapply(x, length))
})



# Find genes with motif in promoter
promoters <- trim(tssGr+promoter_size)
#motifPositions_Cell$Ciliated

genesMotif_Cell <- lapply(motifPositions_Cell, function(motifPositions_Cellgl){
  genesMotif_l <- mclapply(setNames(names(motifPositions_Cellgl), names(motifPositions_Cellgl)), function(TF){
    genesWithMotif <- subsetByOverlaps(promoters, motifPositions_Cellgl[[TF]])
    tibble(TF = TF, target = genesWithMotif$symbol)
  #}, mc.cores = 1)
  }, mc.cores = Ncores)
  
  #print(genesMotif_l)
  bind_rows(genesMotif_l) %>% 
    dplyr::mutate(TF = sub("_.*", "", TF)) %>% 
    dplyr::distinct()
})

genesMotif_Cell
sapply(genesMotif_Cell, nrow)
#------------------------------------------------------------------------------#
#             Filter co-expression modules to make regulons                    #
#------------------------------------------------------------------------------#
# Trim out genes that are not a cis target of the TF

regulons_df_Cell <- lapply(genesMotif_Cell, function(genesMotif){
  inner_join(genesMotif, scenic_corrmodules, by=c("TF", "target")) %>% 
    dplyr::group_by(TF, regulation) %>% 
    dplyr::filter(n() >= min_regulon_size)
})
regulons_df_Cell
lapply(regulons_df_Cell, function(x) length(unique(x$TF)))

saveRDS(regulons_df_Cell, file = p_regulons)

