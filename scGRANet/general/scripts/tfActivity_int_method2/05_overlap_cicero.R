p_regulons       <- "results/integrated/TF_activity_method2_5000/tfRegulons_asDF.RDS"
p_regulonAUC <- "results/integrated/TF_activity_method2_5000/regulonAUC.RDS"
p_archrproj          <- "results/scatac/archr/Save-projcovid6/"
genome               <- "hg19"


p_scenic_adjacencies <- "results/scrna/SCENIC/expr_mat.adjacencies.cor.tsv"
#p_scenic_corrmodules <- "results/integrated/TF_activity_method2/tfModules_asDF.RDS"
p_regulons           <- "results/integrated/TF_activity_method2/tfRegulons_asDF_byCell.RDS"
importance_threshold <- 50/100
promoter_size        <- 5000
min_regulon_size     <- 10
Ncores               <- 20
cellWithPeak         <- 0.1 # Minimum fraction of cells for a given cell type, should have this peak to conserve it



regulons_df_Cell <- readRDS(p_regulons)
regulonAUC <- readRDS(p_regulonAUC)

regulons_df_Cell




#------------------------------------------------------------------------------#

library(tidyverse)
library(ArchR)

#addMotifAnnotations()
# change cutoff to 0.05


#------------------------------------------------------------------------------#
#                   Find genes part of regulons                                #
#------------------------------------------------------------------------------#
if (genome == "hg19") {
  library(EnsDb.Hsapiens.v75)
  genesGr <- genes(EnsDb.Hsapiens.v75)
} else if (genome == "hg38") {
  library(EnsDb.Hsapiens.v86)
  genesGr <- genes(EnsDb.Hsapiens.v86)
}
seqlevelsStyle(genesGr) <- 'UCSC'

# Find genes in regulons
regulon_target_genes <- unique(do.call(c, lapply(regulons_df_Cell, function(x) x$target)))
genesGr <- genesGr[genesGr$symbol %in% regulon_target_genes]

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
#names(motifPositions)

regulon_TFs <- unique(do.call(c, lapply(regulons_df_Cell, function(x) x$TF)))
names(motifPositions)
names()

mapTFs <- tibble(motif_archr = names(motifPositions),
                    TF_archr    = sub("_.*", "", names(motifPositions))) %>% 
  dplyr::mutate(TF = regulon_TFs[match(TF_archr, regulon_TFs)]) %>% 
  dplyr::filter(!is.na(TF)) %>% 
  dplyr::distinct() 

# Remove motifs for which TF not found 
motifPositions <- motifPositions[names(motifPositions) %in% mapTFs$motif_archr]

#------------------------------------------------------------------------------#
#                          Overlap Cicero and regulons                         #
#------------------------------------------------------------------------------#

ciceroToGr <- function(connsdf, column){
  connsdf <- as.data.frame(connsdf)
  x <- as.character(connsdf[,column])
  xdf <- do.call(rbind, strsplit(x, "_"))
  xdf <- as.data.frame(xdf)
  colnames(xdf) <- c("chr", "start", "end")
  xgr <- makeGRangesFromDataFrame(xdf)
  #xgr$cicero <- x
  
  xgr@elementMetadata <- DataFrame(connsdf)
  return(xgr)
}


ciceroOverlap <- function(conns, motifpos, regulon, geneRegion){
  motifIDs <- setNames(names(motifpos), names(motifpos))
  
  # Motif in peak1
  peak1 <- ciceroToGr(conns, "Peak1")
  motifpos_peak1 <- lapply(motifIDs, function(motifID){
    peak_motif <- subsetByOverlaps(peak1, motifpos[[motifID]])
    peak_motif$TF <- motifID
    peak_motif
  })
  # Motif in peak2
  peak2 <- ciceroToGr(conns, "Peak2")
  motifpos_peak2 <- lapply(motifIDs, function(motifID){
    peak_motif <- subsetByOverlaps(peak2, motifpos[[motifID]])
    peak_motif$TF <- motifID
    peak_motif
  })
  
  
  # Overlap of genes in regulon with other peak in cicero conns
  motifdf_tf_peak1 <- lapply(motifIDs, function(motifID){
    # Geat cicero position in peak 1 of TF
    peak2_in1 <- ciceroToGr(motifpos_peak1[[motifID]]@elementMetadata, "Peak2")
    regulonTargets <- regulon$target[regulon$TF %in% motifID]
    
    # Find if gene is in cicero peak2
    ciceromatch <- lapply(regulonTargets, function(regulonTarget){
      targetGr <- geneRegion[geneRegion$symbol %in% regulonTarget]
      targetGr <- subsetByOverlaps(peak2_in1, targetGr)
      if (length(targetGr) >=1) {
        targetGr$target <- regulonTarget
      }
      targetGr
    })
    do.call(rbind, lapply(ciceromatch, function(x) as.data.frame(elementMetadata(x))))
  })
  motifdf_tf_peak1 <- bind_rows(motifdf_tf_peak1) %>% 
    mutate(TF_position = "Peak1")
  
  # Overlap of genes in regulon with other peak in cicero conns
  motifdf_tf_peak2 <- lapply(motifIDs, function(motifID){
    # Geat cicero position in peak 1 of TF
    peak1_in2 <- ciceroToGr(motifpos_peak2[[motifID]]@elementMetadata, "Peak1")
    regulonTargets <- regulon$target[regulon$TF %in% motifID]
    
    # Find if gene is in cicero peak2
    ciceromatch <- lapply(regulonTargets, function(regulonTarget){
      targetGr <- geneRegion[geneRegion$symbol %in% regulonTarget]
      targetGr <- subsetByOverlaps(peak1_in2, targetGr)
      if (length(targetGr) >=1) {
        targetGr$target <- regulonTarget
      }
      targetGr
    })
    do.call(rbind, lapply(ciceromatch, function(x) as.data.frame(elementMetadata(x))))
  })
  motifdf_tf_peak2 <- bind_rows(motifdf_tf_peak2) %>% 
    mutate(TF_position = "Peak2")
  
  
  
  
  return(bind_rows(motifdf_tf_peak1, motifdf_tf_peak2))
}


#names(motifPositions) <- mapTFs$TF
regulons_cnns <- ciceroOverlap(conns, head(motifPositions), regulons_df_Cell$Ciliated, tssGr+5000)


head(conns)



regulons_cnns %>% 
  dplyr::select(Peak1, Peak2, coaccess) %>% 
  distinct() %>% 
  mutate(Regulon = TRUE) %>% 
  bind_rows(conns) %>% 
  ggplot(aes(x = coaccess, color=Regulon)) +
  geom_density()


conns
