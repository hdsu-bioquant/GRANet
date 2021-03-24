p_nmfAUCpySCENIC <- "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/scrna/SCENIC/NMF_AdultLiver/NMF_regulonAUC.RDS"
p_nmfregulonAUC  <- "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_AdultLiver/NMF_regulonAUC.RDS"

p_regulonpySCENIC <- "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/scrna/SCENIC/regulons.csv"
p_regulon <- "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/tfRegulons_asDF.RDS"

#------------------------------------------------------------------------------#
#                               Read Data                                      #
#------------------------------------------------------------------------------#

exprs <- readRDS(p_exprs)
nmfAUCpySCENIC <- readRDS(p_nmfAUCpySCENIC)
nmfregulonAUC  <- readRDS(p_nmfregulonAUC)

hnmfregulonAUC <- HMatrix(nmfregulonAUC, k=8)
hnmfAUCpySCENIC <- HMatrix(nmfAUCpySCENIC, k=8)

as.data.frame(t(hnmfAUCpySCENIC)) %>% 
  rownames_to_column("TF")

wnmfregulonAUC <- WMatrix(nmfregulonAUC, k=11)
wnmfAUCpySCENIC <- WMatrix(nmfAUCpySCENIC, k=11)

xreg <- as.data.frame(wnmfregulonAUC) %>% 
  rownames_to_column("TF") %>% 
  mutate(TFsim = sub(" \\(.*", "", TF))
xsce <- as.data.frame(wnmfAUCpySCENIC) %>% 
  rownames_to_column("TF")%>% 
  mutate(TFsim = sub("\\(.*", "", TF))

xreg_scep <- xreg[!xreg$TFsim %in% xsce$TFsim, ]

table(xreg$TFsim %in% xsce$TFsim)

table(xsce$TFsim %in% xreg$TFsim)

xreg %>% 
  mutate(TFsim = sub(" \\(.*", "", TF))



x %>% 
  dplyr::filter(grepl("zeb", TF, ignore.case=TRUE))


#------------------------------------------------------------------------------#
#         Find co-appearance of cell type and TF in publications               #
#------------------------------------------------------------------------------#

install_github("ropensci/rentrez")

library(rentrez)

res <- entrez_search(db = "pubmed", term = "liver")
res

res <- entrez_search(db = "pubmed", term = "Kuppfer")
res
res$count

res <- entrez_search(db = "pubmed", term = xreg$TFsim[1])
res$count
