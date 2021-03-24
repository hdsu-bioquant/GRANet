options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)


p_regulonAUC         <- as.character(args[1])
p_regulonAUC_devCor  <- as.character(args[2])
p_regulonAUCpySCENIC <- as.character(args[3])
p_RegulonStats       <- as.character(args[4]) 
#------------------------------------------------------------------------------#
# p_regulonAUC         <- "results/integrated/TF_activity/regulonAUC.RDS"
# p_regulonAUC_devCor  <- "results/integrated/TF_activity/regulonAUC_devCorrected.RDS"
# p_regulonAUCpySCENIC <- "results/scrna/SCENIC/auc_mtx.RDS"
# p_RegulonStats       <- "results/scrna/SCENIC/auc_mtx.RDS"
#------------------------------------------------------------------------------#

regulonStats <- function(regulon_path){
  regulon <- readRDS(regulon_path)
  x0 <- basename(regulon_path)
  x1 <- paste0("No. Cell: ", ncol(regulon))
  x2 <- paste0("Regulons: ", nrow(regulon))
  x3 <- paste0("Positive Regulated: ", length(grep("\\(\\+\\)", rownames(regulon))))
  x4 <- paste0("Negative Regulated: ", length(grep("\\(\\-\\)", rownames(regulon))))
  regInfo <- c(x0, x1, x2, x3, x4, "")
  return(regInfo)
}


lines <- c(regulonStats(p_regulonAUC),
           regulonStats(p_regulonAUC_devCor),
           regulonStats(p_regulonAUCpySCENIC))


writeLines(lines, p_RegulonStats)

