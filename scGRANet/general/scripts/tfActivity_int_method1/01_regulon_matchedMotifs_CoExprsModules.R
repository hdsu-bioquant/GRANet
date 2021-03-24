options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

p_scenic_adjacencies <- as.character(args[1])
p_archrproj          <- dirname(as.character(args[2]))
p_scenic_corrmodules <- as.character(args[3])
p_regulons           <- as.character(args[4])
genome               <- as.character(args[5])
importance_threshold <- as.numeric(args[6])/100
promoter_size        <- as.numeric(args[7])
min_regulon_size     <- as.numeric(args[8])
#------------------------------------------------------------------------------#
# p_scenic_adjacencies <- "results/scrna/SCENIC/expr_mat.adjacencies.cor.tsv"
# p_archrproj          <- "results/scatac/archr/Save-projcovid6/"
# p_scenic_corrmodules <- "results/scrna/SCENIC/tfModules_asDF.RDS"
# p_regulons           <- "results/integrated/TF_activity/tfRegulons_asDF.RDS"
# genome               <- "hg19"
# importance_threshold <- 50/100
# promoter_size        <- 5000
# min_regulon_size     <- 10
#------------------------------------------------------------------------------#
#p_archrproj          <- "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_lung/results/scatac/archr/ArchR02_MotifMatch/"
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
saveRDS(scenic_corrmodules, p_scenic_corrmodules)

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

# Find genes in co-expression modules
genesGr <- genesGr[genesGr$symbol %in% unique(scenic_corrmodules$target)]

# get TSS
tssGr <- resize(genesGr, width = 1, fix = "start")
seqlevelsStyle(tssGr) <- 'UCSC'
tssGr

#------------------------------------------------------------------------------#
#           Filter TFs that have no motif in the ATACseq data                  #
#------------------------------------------------------------------------------#
archrproj <- loadArchRProject(path = p_archrproj)

archrproj@peakAnnotation

motifPositions <- getPositions(archrproj)
motifPositions
#names(motifPositions)

scenicTFs <- tibble(motif_archr = names(motifPositions),
                    TF_archr    = sub("_.*", "", names(motifPositions))) %>% 
  dplyr::mutate(TF = scenic_corrmodules$TF[match(TF_archr, scenic_corrmodules$TF)]) %>% 
  dplyr::filter(!is.na(TF)) %>% 
  dplyr::distinct() 


# Remove modules where the TF was not found in the ATACseq data
scenic_corrmodules <- scenic_corrmodules %>% 
  dplyr::filter(TF %in% scenicTFs$TF)

# Remove motifs for TF not found 
motifPositions <- motifPositions[names(motifPositions) %in% scenicTFs$motif_archr]

#------------------------------------------------------------------------------#
#                 Find genes with a motif in its promoter.                     #
#------------------------------------------------------------------------------#
# Find peaks between 5k of tss
# To identify genes that are cis regulated by the TF
# And later remove genes that are downstream target of the TF
promoters <- trim(tssGr+promoter_size)
genesMotif_l <- lapply(setNames(names(motifPositions), names(motifPositions)), function(TF){
  
  genesWithMotif <- subsetByOverlaps(promoters, motifPositions[[TF]])
  tibble(TF = TF, target = genesWithMotif$symbol)
})

genesMotif <- bind_rows(genesMotif_l) %>% 
  dplyr::mutate(TF = sub("_.*", "", TF)) %>% 
  dplyr::distinct()
genesMotif
#------------------------------------------------------------------------------#
#             Filter co-expression modules to make regulons                    #
#------------------------------------------------------------------------------#
# Trim out genes that are not a cis target of the TF
regulons_df <- inner_join(genesMotif, scenic_corrmodules, by=c("TF", "target"))
regulons_df %>% 
  dplyr::group_by(TF, regulation) %>% 
  dplyr::summarise(n = n()) %>% 
  summary()

regulons_df <- regulons_df %>% 
  dplyr::group_by(TF, regulation) %>% 
  dplyr::filter(n() >= min_regulon_size)
regulons_df %>% 
  dplyr::group_by(TF, regulation) %>% 
  dplyr::summarise(n = n()) %>% 
  summary()

saveRDS(regulons_df, p_regulons)


# getPositions()
# peakAnnoEnrichment
# x <- ArchR::getMatrixFromProject(archrproj, useMatrix = "MotifMatrix")
# x@assays@data$deviations
# rownames(x)
# colnames(x)
# 
# ArchR::getAvailableMatrices(archrproj)
# ArchR::getPeakAnnotation(archrproj)
# ArchR::getPeakSet(archrproj)
# 
# getMatches(archrproj)
# # ChromVar motif deviations
# dfVarDev <- getVarDeviations(archrproj, name = "MotifMatrix", plot = FALSE)
# dim(dfVarDev)
# #length(unique(dfVarDev$idx))
# 
