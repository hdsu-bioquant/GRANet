options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

p_regulonAUCpySCENIC <- as.character(args[1])
p_regulonAUCpySCENICRDS <- as.character(args[2])
#------------------------------------------------------------------------------#
# p_regulonAUCpySCENIC <- "results/scrna/SCENIC/auc_mtx.csv"
# p_regulonAUCpySCENICR<- "results/scrna/SCENIC/auc_mtx.RDS"
#------------------------------------------------------------------------------#
library(tidyverse)

regulonAUC_pySCENIC <- read_csv(p_regulonAUCpySCENIC) %>% 
  column_to_rownames("Cell") %>% 
  t()
saveRDS(regulonAUC_pySCENIC, p_regulonAUCpySCENICRDS)
