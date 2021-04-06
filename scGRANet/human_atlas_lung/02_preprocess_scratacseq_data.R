#options(echo=TRUE) # if you want see commands in output file
#args <- commandArgs(TRUE)

# p_scenic_adjacencies <- as.character(args[1])
# p_archrproj          <- dirname(as.character(args[2]))
# p_regulons           <- as.character(args[3])
# genome               <- as.character(args[4])
# importance_threshold <- as.numeric(args[5])/100
# promoter_size        <- as.numeric(args[6])
# min_regulon_size     <- as.numeric(args[7])
# Ncores               <- as.numeric(args[8])
# cellWithPeak         <- as.numeric(args[9])
#------------------------------------------------------------------------------#
# p_scenic_adjacencies <- "results/scrna/SCENIC/expr_mat.adjacencies.cor.tsv"
# p_archrproj          <- "results/scatac/archr/Save-projcovid6/"
# #p_scenic_corrmodules <- "results/integrated/TF_activity_method2/tfModules_asDF.RDS"
# p_regulons           <- "results/integrated/TF_activity_method2/tfRegulons_asDF_byCell.RDS"
genome               <- "hg38"
# importance_threshold <- 50/100
# promoter_size        <- 5000
# min_regulon_size     <- 10
# Ncores               <- 20 
# cellWithPeak         <- 0.1 # Minimum fraction of cells for a given cell type, should have this peak to conserve it
#rawFiles <- c("~/projects/mouse_brain_TF_activity/data/scatac/Lung2_62216_barcoded.bam")
rawFiles <- c("/media/ag-cherrmann/projects/10_charite_covid19/subprojects/human_atlas_TF_activity/data/atac/fragments_processed/lung_SM-A62E9_rep1_fragments.txt",
              "/media/ag-cherrmann/projects/10_charite_covid19/subprojects/human_atlas_TF_activity/data/atac/fragments_processed/lung_SM-A8WNH_rep1_fragments.txt",
              "/media/ag-cherrmann/projects/10_charite_covid19/subprojects/human_atlas_TF_activity/data/atac/fragments_processed/lung_SM-ACCPU_rep1_fragments.txt",
              "/media/ag-cherrmann/projects/10_charite_covid19/subprojects/human_atlas_TF_activity/data/atac/fragments_processed/lung_SM-JF1NZ_rep1_fragments.txt")


archerOutdir <- unique(dirname(rawFiles))
#------------------------------------------------------------------------------#

# Bam file did not have the barcode field by default, thus it was converted into a 
# sam file, and the barcode was added using awk:
# Rsamtools::asSam("~/projects/mouse_brain_TF_activity/data/scatac/Lung2_62216.bam")
# awk '{ if($0 ~ "^@") {print $0} else {split($1,a,":"); gsub(/$/, "\tRG:Z:"a[1]); print} }' Lung2_62216.sam > Lung2_62216_barcoded.sam
# and then converted back to a bam file:
# Rsamtools::asBam("~/projects/mouse_brain_TF_activity/data/scatac/Lung2_62216_barcoded.sam")

# Is is better to use samtools:
# samtools view -h in.bam|awk '{ if($0 ~ "^@") {print $0} else {split($1,a,":"); gsub(/RG:Z:[^\t]*/, "RG:Z:"a[1]); print} }'|samtools view -b -o out.bam


# scanBamHeader(rawFiles)
# myBam <- scanBam(rawFiles)


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
