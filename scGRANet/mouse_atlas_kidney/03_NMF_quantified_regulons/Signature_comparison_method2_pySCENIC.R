
params <- list(p_nmfAUCpySCENIC = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/scrna/SCENIC/NMF_whole/NMF_regulonAUC.RDS",
               p_nmfregulonAUC  = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_regulonAUC.RDS",
               K = 11,
               p_cooccurrance = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_TF_tissue_coocurrance.pdf"
)


# params <- list(p_nmfAUCpySCENIC = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/scrna/SCENIC/NMF_sampled/NMF_regulonAUC.RDS",
#                p_nmfregulonAUC  = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_sampled/NMF_regulonAUC.RDS",
#                K = 11
# )
# p_regulonpySCENIC <- "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/scrna/SCENIC/regulons.csv"
# p_regulon <- "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/tfRegulons_asDF.RDS"

#------------------------------------------------------------------------------#
#                               Read Data                                      #
#------------------------------------------------------------------------------#

#exprs <- readRDS(p_exprs)
nmfregulonAUC  <- readRDS(params$p_nmfregulonAUC)
nmfAUCpySCENIC <- readRDS(params$p_nmfAUCpySCENIC)

hnmfregulonAUC <- HMatrix(nmfregulonAUC, k=params$K)
hnmfAUCpySCENIC <- HMatrix(nmfAUCpySCENIC, k=params$K)

hreg <- as.data.frame(hnmfregulonAUC) %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("CellID")
hsce <- as.data.frame(hnmfAUCpySCENIC) %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("CellID")


wnmfregulonAUC <- WMatrix(nmfregulonAUC, k=params$K)
wnmfAUCpySCENIC <- WMatrix(nmfAUCpySCENIC, k=params$K)

wreg <- as.data.frame(wnmfregulonAUC) %>% 
  rownames_to_column("TF") %>% 
  mutate(TFsim = sub(" \\(.*", "", TF))
wsce <- as.data.frame(wnmfAUCpySCENIC) %>% 
  rownames_to_column("TF")%>% 
  mutate(TFsim = sub("\\(.*", "", TF))

# table(xreg$TFsim %in% xsce$TFsim)
# table(xsce$TFsim %in% xreg$TFsim)
#xreg_scep <- xreg[!xreg$TFsim %in% xsce$TFsim, ]

#------------------------------------------------------------------------------#
#         Find co-appearance of cell type and TF in publications               #
#------------------------------------------------------------------------------#
#install_github("ropensci/rentrez")
library(rentrez)

# res_tissue <- entrez_search(db = "pubmed", term = "brain")
# res_tissue

pubmed_signature <- function(w_matrix, signatureIndex, tissueName) {
  signa <-  w_matrix[,signatureIndex]
  signa <-  sort(signa, decreasing = TRUE)
  TFnames <-  sub(" \\(.*|\\(.*","",names(signa))
  searchTerm <- paste0("(", TFnames, ") AND (", tissueName, ")")
  names(searchTerm) <- names(signa)
  pubmed_search_list <- lapply(searchTerm, function(t) entrez_search(db = "pubmed", term = t))
  return(pubmed_search_list)
}

cooccurrence_gg <- function(wPubmedregulonAUC, wPubmedpySCENIC, tissue){
  trank <- 0.2
  tcount <- 0.90
  
  dfm2 <- tibble(RegulonID = names(wPubmedregulonAUC),
                 rank = 1:length(wPubmedregulonAUC),
                 count = sapply(wPubmedregulonAUC, function(x) x$count)) %>%
    mutate(TFnames = sub(" \\(.*|\\(.*","", RegulonID)) %>% 
    dplyr::filter(!grepl("\\(-\\)", RegulonID)) %>% 
    mutate(count = log(count)) 
  
  dfps <- tibble(RegulonID = names(wPubmedpySCENIC),
                 rank = 1:length(wPubmedpySCENIC),
                 count = sapply(wPubmedpySCENIC, function(x) x$count)) %>% 
    mutate(TFnames = sub(" \\(.*|\\(.*","", RegulonID)) %>% 
    mutate(count = log(count))
  
  a <- dfm2 %>% 
    mutate(Specific = !TFnames %in% dfps$TFnames) %>% 
    ggplot(aes(x = rank, y = count)) +
    geom_point(aes(color = Specific)) +
    ggrepel::geom_text_repel(data = dfm2 %>% filter(rank < quantile(rank, trank) & count > quantile(count, tcount)),
              aes(label=RegulonID,hjust=0,vjust=0)) +
    ggtitle("Integrative inference of regulon activity") +
    xlab(paste0("Regulon Rank\nNMF ", tissue, " signature")) +
    ylab("Co-occurrence of TF and tissue in indexed publications") +
    cowplot::theme_cowplot()
  
  b <- dfps %>% 
    mutate(Specific = !TFnames %in% dfm2$TFnames) %>% 
    ggplot(aes(x = rank, y = count)) +
    geom_point(aes(color = Specific)) +
    ggrepel::geom_text_repel(data = dfps %>% filter(rank < quantile(rank, trank) & count > quantile(count, tcount)),
                             aes(label=RegulonID,hjust=0,vjust=0)) +
    ggtitle("pySCENIC") +
    xlab(paste0("Regulon Rank\nNMF ", tissue, " signature")) +
    ylab("Co-occurrence of TF and tissue in indexed publications") +
    cowplot::theme_cowplot()
  
  return( a + b )
}

SSFintRegul <- SignatureSpecificFeatures(nmfregulonAUC, k=params$K)
SSFpySCENIC <- SignatureSpecificFeatures(nmfAUCpySCENIC, k=params$K)

SSFintRegul <- SSFintRegul[c(6,10,9,2)]
SSFpySCENIC <- SSFpySCENIC[c(1,3,10,5)]
names(SSFintRegul) <- names(SSFpySCENIC) <- c("Brain", "Thymus", "Kidney", "Liver")



cooccurrence_gg2 <- function(wPubmedregulonAUC, wPubmedpySCENIC, tissue){
  trank <- 0.2
  tcount <- 0.90
  
  dfm2 <- tibble(RegulonID = names(wPubmedregulonAUC),
                 rank = 1:length(wPubmedregulonAUC),
                 count = sapply(wPubmedregulonAUC, function(x) x$count)) %>%
    mutate(TFnames = sub(" \\(.*|\\(.*","", RegulonID)) %>% 
    dplyr::filter(!grepl("\\(-\\)", RegulonID)) %>% 
    mutate(count = log(count)) 
  
  dfps <- tibble(RegulonID = names(wPubmedpySCENIC),
                 rank = 1:length(wPubmedpySCENIC),
                 count = sapply(wPubmedpySCENIC, function(x) x$count)) %>% 
    mutate(TFnames = sub(" \\(.*|\\(.*","", RegulonID)) %>% 
    mutate(count = log(count))
  
  #print(dfm2)
  
  a <- dfm2 %>% 
    mutate(Specific = !TFnames %in% dfps$TFnames) %>% 
    ggplot(aes(x = rank, y = count)) +
    geom_point(aes(color = Specific)) +
    
    ggrepel::geom_text_repel(data = dfm2 %>% dplyr::filter(RegulonID %in% SSFintRegul[[tissue]]),
                             aes(label=RegulonID,hjust=0,vjust=0)) +
    ggtitle("Integrative inference of regulon activity") +
    xlab(paste0("Regulon Rank\nNMF ", tissue, " signature")) +
    ylab("Co-occurrence of TF and tissue in indexed publications") +
    cowplot::theme_cowplot()
  
  b <- dfps %>% 
    mutate(Specific = !TFnames %in% dfm2$TFnames) %>% 
    ggplot(aes(x = rank, y = count)) +
    geom_point(aes(color = Specific)) +
    ggrepel::geom_text_repel(data = dfps %>% dplyr::filter(RegulonID %in% SSFpySCENIC[[tissue]]),
                             aes(label=RegulonID,hjust=0,vjust=0)) +
    ggtitle("pySCENIC") +
    xlab(paste0("Regulon Rank\nNMF ", tissue, " signature")) +
    ylab("Co-occurrence of TF and tissue in indexed publications") +
    cowplot::theme_cowplot()
  
  return( a + b )
}
#cooccurrence_gg2(wPubmedregulonAUC, wPubmedpySCENIC, "Brain")


cooccurrence_plist <- list()
# Search Brain signature
# Signature 6 for method 2
# Signature 1 for pySCENIC
wPubmedregulonAUC <- pubmed_signature(wnmfregulonAUC, signatureIndex=6, tissueName = "brain") 
wPubmedpySCENIC   <- pubmed_signature(wnmfAUCpySCENIC, signatureIndex=1, tissueName = "brain") 
cooccurrence_plist[[1]] <- cooccurrence_gg2(wPubmedregulonAUC, wPubmedpySCENIC, "Brain")
cooccurrence_plist[[1]]

# Search Thymus signature
# Signature 10 for method 2
# Signature 3 for pySCENIC
thymuswPubmedregulonAUC <- pubmed_signature(wnmfregulonAUC, signatureIndex=10, tissueName = "thymus") 
thymuswPubmedpySCENIC   <- pubmed_signature(wnmfAUCpySCENIC, signatureIndex=3, tissueName = "thymus") 
cooccurrence_plist[[2]] <- cooccurrence_gg2(thymuswPubmedregulonAUC, thymuswPubmedpySCENIC, "Thymus")
cooccurrence_plist[[2]]

# Search kidney signature
# Signature 9 for method 2
# Signature 10 for pySCENIC
kidneywPubmedregulonAUC <- pubmed_signature(wnmfregulonAUC, signatureIndex=9, tissueName = "kidney") 
kidneywPubmedpySCENIC   <- pubmed_signature(wnmfAUCpySCENIC, signatureIndex=10, tissueName = "kidney") 
cooccurrence_plist[[3]] <- cooccurrence_gg2(kidneywPubmedregulonAUC, kidneywPubmedpySCENIC, "Kidney")
cooccurrence_plist[[3]]

# Search liver signature
# Signature 2 for method 2
# Signature 5 for pySCENIC
liverwPubmedregulonAUC <- pubmed_signature(wnmfregulonAUC, signatureIndex=2, tissueName = "liver") 
liverwPubmedpySCENIC   <- pubmed_signature(wnmfAUCpySCENIC, signatureIndex=5, tissueName = "liver") 
cooccurrence_plist[[4]] <- cooccurrence_gg2(liverwPubmedregulonAUC, liverwPubmedpySCENIC, "Liver")
cooccurrence_plist[[4]]


pdf(file = params$p_cooccurrance, width=15, height=7)
for (p in cooccurrence_plist) {
  plot(p)
}
dev.off()





#------------------------------------------------------------------------------#
#                  Wordcloud abstract top regulons                             #
#------------------------------------------------------------------------------#
pubmed_signature_top <- function(w_matrix, signatureIndex, tissueName, topn) {
  signa <-  w_matrix[,signatureIndex]
  signa <-  sort(signa, decreasing = TRUE)
  TFnames <-  sub(" \\(.*|\\(.*","",names(signa))
  searchTerm <- paste0("(", TFnames, ") AND (", tissueName, ")")
  #searchTerm <- paste0("(", TFnames, ")")
  names(searchTerm) <- names(signa)
  searchTerm <- searchTerm[1:topn]
  pubmed_search_list <- lapply(searchTerm, function(t) entrez_search(db = "pubmed", term = t, use_history=TRUE))
  return(pubmed_search_list)
}

pubmed_parse_abstract <- function(pubmedlist) {
  abslist <- lapply(pubmedlist, function(x){
    print(x)
    recs <- entrez_fetch(db="pubmed", web_history=x$web_history, retmax=1000000, rettype="xml", parsed=TRUE)
    recs <- parse_pubmed_xml(recs)
    # list(doi      = unlist(lapply(recs, "[[", "doi")),
    #      abstract = unlist(lapply(recs, function(y) paste(y$abstract, collapse = " "))))
    list(doi = unlist(lapply(recs, function(x) {
      d <- x$doi
      ifelse(is_character(d), d, NA)
      })),
      abstract = lapply(recs, function(y) paste(y$abstract, collapse = " "))
    )
  })
  #return(abslist)
  abslist <- list(doi = do.call(c, lapply(abslist, "[[", "doi")),
                  abstract = do.call(c, lapply(abslist, "[[", "abstract")))
  
  # Remove duplicated papers
  mask <- duplicated(abslist$doi) | is.na(abslist$doi)
  #return(abslist$abstract[!mask])
  
  x <- abslist$abstract[!mask]
  # Create a corpus with abstract
  docs <- Corpus(VectorSource(x))
  # Clean the text data
  docs <- docs %>%
    tm_map(removeNumbers) %>%
    tm_map(removePunctuation) %>%
    tm_map(stripWhitespace)
  docs <- tm_map(docs, content_transformer(tolower))
  docs <- tm_map(docs, removeWords, stopwords("english"))
  # Create a document-term-matrix 
  dtm <- TermDocumentMatrix(docs) 
  word_matrix <- as.matrix(dtm) 
  words <- sort(rowSums(word_matrix),decreasing=TRUE) 
  return(data.frame(word = names(words),freq=words) )
}

# Search Brain signature
# Signature 6 for method 2
psBrain_intRegul <- pubmed_signature_top(wnmfregulonAUC, 6, tissueName="brain", topn=10)
wordsBrain_intRegul <- pubmed_parse_abstract(psBrain_intRegul)

# Search Thymus signature
# Signature 10 for method 2
psthymus_intRegul <- pubmed_signature_top(wnmfregulonAUC, 10, tissueName="thymus", topn=10)
wordsThymus_intRegul <- pubmed_parse_abstract(psthymus_intRegul)

# Search kidney signature
# Signature 9 for method 2
pskidney_intRegul <- pubmed_signature_top(wnmfregulonAUC, 9, tissueName="kidney", topn=10)
wordsKidney_intRegul <- pubmed_parse_abstract(pskidney_intRegul)

# Search liver signature
# Signature 2 for method 2
psliver_intRegul <- pubmed_signature_top(wnmfregulonAUC, 2, tissueName="liver", topn=10)
wordsLiver_intRegul <- pubmed_parse_abstract(psliver_intRegul)

# Search bone marrow signature
# Signature 8 for method 2
psBoneMarrow_intRegul <- pubmed_signature_top(wnmfregulonAUC, 8, tissueName="bone marrow", topn=10)
wordsBoneMarrow_intRegul <- pubmed_parse_abstract(psBoneMarrow_intRegul)

# Search Lung signature
# Signature 1 for method 2
psLung_intRegul <- pubmed_signature_top(wnmfregulonAUC, 1, tissueName="Lung", topn=10)
wordsLung_intRegul <- pubmed_parse_abstract(psLung_intRegul)

# Search SmallIntestine signature
# Signature 3 for method 2
psSmallIntestine_intRegul <- pubmed_signature_top(wnmfregulonAUC, 3, tissueName="small intestine", topn=10)
wordsSmallIntestine_intRegul <- pubmed_parse_abstract(psSmallIntestine_intRegul)

# Search Spleen signature
# Signature 7 for method 2
psSpleen_intRegul <- pubmed_signature_top(wnmfregulonAUC, 7, tissueName="spleen", topn=10)
wordsSpleen_intRegul <- pubmed_parse_abstract(psSpleen_intRegul)

# Search Testis signature
# Signature 4 for method 2
psTestis_intRegul <- pubmed_signature_top(wnmfregulonAUC, 4, tissueName="testis", topn=10)
wordsTestis_intRegul <- pubmed_parse_abstract(psTestis_intRegul)

words_df_list <- list(
  Brain = wordsBrain_intRegul,
  Thymus = wordsThymus_intRegul,
  Kidney = wordsKidney_intRegul,
  Liver = wordsLiver_intRegul,
  BoneMarrow = wordsBoneMarrow_intRegul,
  Lung = wordsLung_intRegul,
  SmallIntestine = wordsSmallIntestine_intRegul,
  Spleen = wordsSpleen_intRegul,
  Testis = wordsTestis_intRegul
  )

# save it in html
library(htmlwidgets)
library(webshot)
# webshot::install_phantomjs()


lapply(names(words_df_list), function(id){
  words_df <- words_df_list[[id]]
  df <- words_df %>% dplyr::filter(!word %in% c("cell", "cells", "gene", "genes", "expression"))
  my_graph <- wordcloud2(data=df, shape="circle")
  saveWidget(my_graph, "tmp/tmp.html", selfcontained = FALSE)
  # and in png or pdf
  p <- file.path(dirname(params$p_nmfregulonAUC), paste0("wordCloud_", id, ".pdf"))
  size = 1500
  webshot("tmp/tmp.html",p, delay =5, vwidth = size, vheight=size)
})



# lapply(words_df_list, function(words_df){
#   df <- words_df %>% dplyr::filter(!word %in% c("cell", "cells", "gene", "genes"))
#   my_graph <- wordcloud2(data=df, shape="circle")
#   saveWidget(my_graph, "tmp/tmp.html", selfcontained = FALSE)
#   # and in png or pdf
#   webshot("tmp/tmp.html","fig_1.pdf", delay =5, vwidth = 480, vheight=480)
# })

getwd()

# and in png or pdf
webshot("tmp.html","fig_1.pdf", delay =5, vwidth = 480, vheight=480)


# 
# df <- wordsThymus_intRegul
# 
# 
# 
# 
# # Generate the word cloud
# wordcloud(words = df$word, freq = df$freq, min.freq = 1, 
#           max.words=200, random.order=FALSE, rot.per=0.35,
#           colors=brewer.pal(8, "Dark2"))
# 
# wordcloud2(data=df, shape="circle")
# 
# 
# 
# 
# 
# 
# unlist(lapply(x$`Sox10 (+)`$doi, function(d) ifelse(is_character(d), d, NA)))
# 
# 
# 
# abstractTF_intRegul <- mclapply(pubmedTF_intRegul, function(x){
#   recs <- entrez_fetch(db="pubmed", web_history=x$web_history, retmax=1000000, rettype="xml", parsed=TRUE)
#   recs <- parse_pubmed_xml(recs)
#   #print(length(recs))
#   #print(unlist(lapply(recs, "[[", "doi")))
#   #print(lapply(recs, "[[", "abstract"))
#   list(doi      = unlist(lapply(recs, "[[", "doi")),
#        abstract = unlist(lapply(recs, function(y) paste(y$abstract, collapse = " "))))
# }, mc.cores = 3)
# 
# 
# 
# wPubmedpySCENIC
# 
# 
# 
# 
# 
# pubmedTF_intRegul <- lapply(setNames(sub(" \\(.*|\\(.*","", rownames(wnmfregulonAUC)), 
#                                      rownames(wnmfregulonAUC)), 
#                             function(t) entrez_search(db="pubmed", term=t, use_history=TRUE))
# pubmedTF_pySCENIC <- lapply(setNames(sub(" \\(.*|\\(.*","", rownames(wnmfAUCpySCENIC)), 
#                                      rownames(wnmfAUCpySCENIC)), 
#                             function(t) entrez_search(db="pubmed", term=t, use_history=TRUE))
# 
# abstractTF_intRegul <- mclapply(pubmedTF_intRegul, function(x){
#   recs <- entrez_fetch(db="pubmed", web_history=x$web_history, retmax=1000000, rettype="xml", parsed=TRUE)
#   recs <- parse_pubmed_xml(recs)
#   #print(length(recs))
#   #print(unlist(lapply(recs, "[[", "doi")))
#   #print(lapply(recs, "[[", "abstract"))
#   list(doi      = unlist(lapply(recs, "[[", "doi")),
#          abstract = unlist(lapply(recs, function(y) paste(y$abstract, collapse = " "))))
# }, mc.cores = 3)
# 
# 
# 
# 
# table(sapply(abstractTF_intRegul, function(x) class(x)))
# 
# 
# abstractTF_intRegul
# 
# class(abstractTF_intRegul[["Zbtb7a (-)"]])
# 
# recs <- entrez_fetch(db="pubmed", web_history=pubmedTF_intRegul[[1]]$web_history, retmax=1000000, rettype="xml", parsed=TRUE)
# parse_pubmed_xml(recs)
# class(recs)
# rec[[1]]$doi
# 
# 
# pubmedids <- unique(unlist(lapply(wPubmedregulonAUC[1:20], "[[", "ids")))
# #entrez_fetch(pubmedids, db="pubmed", parsed = TRUE, rettype = "xml")
# rec <- parse_pubmed_xml(entrez_fetch(pubmedids, db="pubmed", parsed = TRUE, rettype = "xml"))
# rec[[1]]$abstract
# 
# unlist(lapply(rec, "[[", "abstract"))
# 
# rec$doi
# 
# 
# 
# 
# 
# 
# # Create a corpus with abstract
# docs <- Corpus(VectorSource(unlist(lapply(rec, "[[", "abstract"))))
# # Clean the text data
# docs <- docs %>%
#   tm_map(removeNumbers) %>%
#   tm_map(removePunctuation) %>%
#   tm_map(stripWhitespace)
# docs <- tm_map(docs, content_transformer(tolower))
# docs <- tm_map(docs, removeWords, stopwords("english"))
# # Create a document-term-matrix 
# dtm <- TermDocumentMatrix(docs) 
# word_matrix <- as.matrix(dtm) 
# words <- sort(rowSums(word_matrix),decreasing=TRUE) 
# df <- data.frame(word = names(words),freq=words)
# # Generate the word cloud
# wordcloud(words = df$word, freq = df$freq, min.freq = 1, 
#           max.words=200, random.order=FALSE, rot.per=0.35,
#           colors=brewer.pal(8, "Dark2"))
# 
# wordcloud2(data=df, shape="circle")
# 
# 
# 
# #install.packages("wordcloud")
# install.packages("wordcloud2")
# library(wordcloud2)
# library(tm)
# 
# Corpus(VectorSource(text))
# 
# 
# a <- tibble(RegulonID = names(wPubmedregulonAUC),
#        rank = 1:length(wPubmedregulonAUC),
#        count = sapply(wPubmedregulonAUC, function(x) x$count)) %>%
#   dplyr::filter(!grepl("\\(-\\)", RegulonID)) %>% 
#   mutate(count = ifelse(count > quantile(count, 0.9), as.numeric(quantile(count, 0.9)), count)) %>% 
#   ggplot(aes(x = rank, y = count)) +
#   geom_point() +
#   cowplot::theme_cowplot()
# 
# dfm2 <- tibble(RegulonID = names(wPubmedregulonAUC),
#        rank = 1:length(wPubmedregulonAUC),
#        count = sapply(wPubmedregulonAUC, function(x) x$count)) %>%
#   mutate(TFnames = sub(" \\(.*|\\(.*","", RegulonID)) %>% 
#   dplyr::filter(!grepl("\\(-\\)", RegulonID)) %>% 
#   mutate(count = log(count)) 
# 
# dfps <- tibble(RegulonID = names(wPubmedpySCENIC),
#        rank = 1:length(wPubmedpySCENIC),
#        count = sapply(wPubmedpySCENIC, function(x) x$count)) %>% 
#   mutate(TFnames = sub(" \\(.*|\\(.*","", RegulonID)) %>% 
#   mutate(count = log(count))
# 
# dfm2 %>% 
#   mutate(Specific = !TFnames %in% dfps$TFnames) %>% 
#   ggplot(aes(x = rank, y = count, color = Specific)) +
#   geom_point() +
#   cowplot::theme_cowplot()
# 
# 
# dfps %>% 
#   mutate(Specific = !TFnames %in% dfm2$TFnames) %>% 
#   ggplot(aes(x = rank, y = count, color = Specific)) +
#   geom_point() +
#   cowplot::theme_cowplot()
# 
# 
# 
# tibble(RegulonID = names(wPubmedpySCENIC),
#        rank = 1:length(wPubmedpySCENIC),
#        count = sapply(wPubmedpySCENIC, function(x) x$count)) %>% 
#   mutate(count = log(count)) %>% 
#   ggplot(aes(x = rank, y = count)) +
#   geom_point() +
#   cowplot::theme_cowplot()
# 
# b <- tibble(RegulonID = names(wPubmedpySCENIC),
#        rank = 1:length(wPubmedpySCENIC),
#        count = sapply(wPubmedpySCENIC, function(x) x$count)) %>% 
#   mutate(count = ifelse(count > quantile(count, 0.9), as.numeric(quantile(count, 0.9)), count)) %>% 
#   ggplot(aes(x = rank, y = count)) +
#   geom_point() +
#   cowplot::theme_cowplot()
# 
# a + b
# 
# as.numeric(quantile(1:100, 0.9))
# 
# 
# plot(sapply(wPubmedregulonAUC, function(tfpubmed){
#   sum(tfpubmed$ids %in% res_tissue$ids)
# }))
# 
# 
# res <- entrez_search(db = "pubmed", term = "Kuppfer")
# res
# res$count
# 
# res <- entrez_search(db = "pubmed", term = xreg$TFsim[1])
# res$count
