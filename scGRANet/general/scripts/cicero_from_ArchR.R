library(cicero)

# p_scenic_adjacencies <- "results/scrna/SCENIC/expr_mat.adjacencies.cor.tsv"
# p_scenic_corrmodules <- "results/scrna/SCENIC/tfModules_asDF.RDS"
# p_regulons           <- "results/integrated/TF_activity/tfRegulons_asDF.RDS"
# genome               <- "hg19"
# importance_threshold <- 50/100
# promoter_size        <- 5000
# min_regulon_size     <- 10

p_archrproj <- "results/scatac/archr/Save-projcovid6/"
p_cicerores <- "results/scatac/cicero.RDS"

archrproj <- loadArchRProject(path = p_archrproj)
ArchR::getAvailableMatrices(archrproj)


peakMatrix <- ArchR::getMatrixFromProject(archrproj, useMatrix = "PeakMatrix")
dim(peakMatrix)
peakMatrix@assays@data$PeakMatrix
peakMatrix@colData
rownames(x)
colnames(x)



#------------------------------------------------------------------------------#
#                           Cicero from ArchR                                  #
#------------------------------------------------------------------------------#

# peaks matrix data 
indata <- peakMatrix@assays@data$PeakMatrix
# format cell info
cellinfo <- as.data.frame(peakMatrix@colData)
# format peak info
peaks <- ArchR::getPeakSet(archrproj)
peaks$Cluster <- names(peaks)
names(peaks) <- NULL
peakinfo <- as.data.frame(peaks) %>% 
  mutate(site_name = paste(seqnames, start, end, sep="_")) %>% 
  column_to_rownames("site_name") %>% 
  as.data.frame()

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

# make CDS
fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
pd <- methods::new("AnnotatedDataFrame", data = cellinfo)
input_cds <-  suppressWarnings(newCellDataSet(indata,
                                              phenoData = pd,
                                              featureData = fd,
                                              expressionFamily=VGAM::binomialff(),
                                              lowerDetectionLimit=0))
input_cds@expressionFamily@vfamily <- "binomialff"
input_cds <- monocle::detectGenes(input_cds)

#Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 


#------------------------------------------------------------------------------#
#                           Create a Cicero CDS                                #
#------------------------------------------------------------------------------#

set.seed(2017)
input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)

# *** if you are using Monocle 3 alpha, you need to run the following line as well!
# NOTE: Cicero does not yet support the Monocle 3 beta release (monocle3 package). We hope
# to update soon!
#input_cds <- preprocessCDS(input_cds, norm_method = "none")
input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                             reduction_method = 'tSNE', norm_method = "none")




# Next, we access the tSNE coordinates from the input CDS object where they are 
# stored by Monocle and run make_cicero_cds:
tsne_coords <- t(reducedDimA(input_cds))
row.names(tsne_coords) <- row.names(pData(input_cds))
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)


#------------------------------------------------------------------------------#
#                                 Run Cicero                                   #
#------------------------------------------------------------------------------#

data("human.hg19.genome")
sample_genome <- subset(human.hg19.genome, V1 == "chr18")
conns <- run_cicero(cicero_cds, sample_genome) # Takes a few minutes to run
head(conns)

saveRDS(list(input_cds =input_cds,
             cicero_cds = cicero_cds,
             conns = conns), p_cicerores)



# compare regulons to cicero connection
#get gene position and subsetoverlap conns
# get motif position and subsetoverlap
