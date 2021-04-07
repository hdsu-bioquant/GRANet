
params <- list(p_nmfAUCpySCENIC = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/scrna/SCENIC/NMF_whole/NMF_regulonAUC.RDS",
               p_nmfregulonAUC  = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_regulonAUC.RDS",
               K = 11,
               p_cooccurrance = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_TF_tissue_coocurrance.pdf",
               p_exprs = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/data/scrna/Han_atlas/MCA_Adult_9Tissues_Seurat.RDS"
)

#------------------------------------------------------------------------------#
#                               Read Data                                      #
#------------------------------------------------------------------------------#
library(ButchR)
library(wordcloud)

nmfregulonAUC  <- readRDS(params$p_nmfregulonAUC)

#------------------------------------------------------------------------------#
#                     Signature specific regulons                              #
#------------------------------------------------------------------------------#
SSFmintRegul <- SignatureSpecificFeatures(nmfregulonAUC, k=params$K, return_all_features=TRUE)
SSFmintRegul <- SSFmintRegul[,c(6,10,9,2,8,1,3,7,4)]
colnames(SSFmintRegul) <- c("Brain", "Thymus", "Kidney", "Liver", "BoneMarrow",
                            "Lung", "SmallIntestine", "Spleen", "Testis")
SSFmintRegul <- SSFmintRegul[rowSums(SSFmintRegul) == 1,]

#------------------------------------------------------------------------------#
#                            W matrix.                                         #
#------------------------------------------------------------------------------#

WnmfregulonAUC <- WMatrix(nmfregulonAUC, k=params$K)
WnmfregulonAUC <- WnmfregulonAUC[,c(6,10,9,2,8,1,3,7,4)]
colnames(WnmfregulonAUC) <- c("Brain", "Thymus", "Kidney", "Liver", "BoneMarrow",
                            "Lung", "SmallIntestine", "Spleen", "Testis")

WnmfregulonAUC <- WnmfregulonAUC[rownames(SSFmintRegul),]

Heatmap(WnmfregulonAUC/rowMaxs(WnmfregulonAUC), col = inferno(100), show_column_names = TRUE, show_row_names = FALSE)


apply(SSFmintRegul, 2, function(x){
  x <- data.frame(TF = names(x[x == 1])) %>% 
    dplyr::filter(!grepl(" \\(-.*", TF)) %>% 
    mutate(W =  rowMaxs(WnmfregulonAUC[match(TF, rownames(WnmfregulonAUC)),])) %>% 
    mutate(TF = sub(" \\(.*", "", TF )) %>% 
    arrange(-W) %>% 
    mutate(f = rev(1:n()))
  
  wordcloud::wordcloud(words = x$TF, freq = x$f, min.freq = 1, 
            max.words=200, random.order=FALSE, rot.per=0.35,
            colors=brewer.pal(8, "Dark2"))
})




pdf(file = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_WMatrix.pdf"
    , width=10, height=7)

Heatmap(WnmfregulonAUC/rowMaxs(WnmfregulonAUC), col = inferno(100), show_column_names = TRUE, show_row_names = FALSE)

Heatmap(WnmfregulonAUC/rowMaxs(WnmfregulonAUC), col = inferno(100), show_column_names = TRUE, show_row_names = FALSE, cluster_columns = FALSE)

apply(SSFmintRegul, 2, function(x){
  x <- data.frame(TF = names(x[x == 1])) %>% 
    dplyr::filter(!grepl(" \\(-.*", TF)) %>% 
    mutate(W =  rowMaxs(WnmfregulonAUC[match(TF, rownames(WnmfregulonAUC)),])) %>% 
    mutate(TF = sub(" \\(.*", "", TF )) %>% 
    arrange(-W) %>% 
    mutate(f = rev(1:n()))
  
  wordcloud::wordcloud(words = x$TF, freq = x$f, min.freq = 1, 
                       max.words=200, random.order=FALSE, rot.per=0.35,
                       colors=brewer.pal(8, "Dark2"))
})


dev.off()





