
params <- list(p_nmfAUCpySCENIC = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/scrna/SCENIC/NMF_whole/NMF_regulonAUC.RDS",
               p_nmfregulonAUC  = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_regulonAUC.RDS",
               K = 11,
               p_cooccurrance = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/results/integrated/TF_activity_method2_500000/NMF_whole/NMF_TF_tissue_coocurrance.pdf"
)


#------------------------------------------------------------------------------#
#                               Read Data                                      #
#------------------------------------------------------------------------------#

#exprs <- readRDS(p_exprs)
nmfregulonAUC  <- readRDS(params$p_nmfregulonAUC)
nmfAUCpySCENIC <- readRDS(params$p_nmfAUCpySCENIC)


#------------------------------------------------------------------------------#
#                     Specific regulon for each tissue                         #
#------------------------------------------------------------------------------#

# Search Brain signature          # Signature 6 for method 2
# Search Thymus signature         # Signature 10 for method 2
# Search kidney signature         # Signature 9 for method 2
# Search liver signature          # Signature 2 for method 2
# Search bone marrow signature    # Signature 8 for method 2
# Search Lung signature           # Signature 1 for method 2
# Search SmallIntestine signature # Signature 3 for method 2
# Search Spleen signature         # Signature 7 for method 2
# Search Testis signature         # Signature 4 for method 2

SSFintRegul <- SignatureSpecificFeatures(nmfregulonAUC, k=params$K, return_all_features=FALSE)
SSFintRegul <- SSFintRegul[c(6,10,9,2,8,1,3,7,4)]
names(SSFintRegul) <- c("Brain", "Thymus", "Kidney", "Liver", "BoneMarrow",
                            "Lung", "SmallIntestine", "Spleen", "Testis")
# Get only name
SSFintRegul <- lapply(SSFintRegul, function(x){
  x <- x[!grepl("\\(-\\)", x)]
  sub(" \\(.*", "", x)
})
SSFintRegul

#------------------------------------------------------------------------------#
#         Find co-appearance of cell type and TF in publications               #
#------------------------------------------------------------------------------#
#install_github("ropensci/rentrez")
library(rentrez)


#------------------------------------------------------------------------------#
#                  Wordcloud abstract SSF regulons                             #
#------------------------------------------------------------------------------#
pubmed_ssf <- function(TFnames, tissueName) {
  searchTerm <- paste0("(", TFnames, ") AND (", tissueName, ")")
  #searchTerm <- paste0("(", TFnames, ")")
  #searchTerm <- paste0('(', TFnames, ') AND "has abstract"[Filter]' )
  names(searchTerm) <- TFnames
  pubmed_search_list <- lapply(searchTerm, function(t) entrez_search(db = "pubmed", term = t, use_history=TRUE))
  return(pubmed_search_list)
}

pubmed_parse_abstract <- function(pubmedlist) {
  abslist <- lapply(pubmedlist, function(x){
    print(x)
    recs <- entrez_fetch(db="pubmed", web_history=x$web_history, retmax=1000000, rettype="xml", parsed=TRUE)
    recs <- parse_pubmed_xml(recs)
    #print(recs)
    # list(doi      = unlist(lapply(recs, "[[", "doi")),
    #      abstract = unlist(lapply(recs, function(y) paste(y$abstract, collapse = " "))))
    list(doi = unlist(lapply(recs, function(x) {
      if (is.character(x)) {
        #print(x)
        d <- NA
      } else {
        d <- x$doi
      }
      ifelse(is_character(d), d, NA)
    })),
    abstract = lapply(recs, function(y) {
      if (is.character(y)) {
        NA
      } else {
        paste(y$abstract, collapse = " ")
      }
      
    })
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


tissueIDs <- setNames(c("Brain", "Thymus", "Kidney", "Liver", "Bone Marrow",
                        "Lung", "Small Intestine", "Spleen", "Testis"),
                      c("Brain", "Thymus", "Kidney", "Liver", "BoneMarrow",
                        "Lung", "SmallIntestine", "Spleen", "Testis"))

#pubmedSearch_SSFintRegul <- lapply(SSFintRegul, pubmed_ssf)
pubmedSearch_SSFintRegul <- lapply(setNames(names(SSFintRegul), names(SSFintRegul)), function(tissue){
  pubmed_ssf(SSFintRegul[[tissue]], tissueIDs[tissue])
})


# Remove 0 count searchs
pubmedSearch_SSFintRegul <- lapply(pubmedSearch_SSFintRegul, function(x) {
  #sum(!sapply(x, "[[", "count") > 0)
  idx <- sapply(x, "[[", "count") > 0
  x[idx]
})
sapply(pubmedSearch_SSFintRegul, function(x) {
  sum(!sapply(x, "[[", "count") > 0)
})

words_SSFintRegul <- lapply(pubmedSearch_SSFintRegul, pubmed_parse_abstract)


# # Search Brain signature
# wordsBrain_intRegul <- pubmed_parse_abstract(pubmedSearch_SSFintRegul$Brain)
# # Search Thymus signature
# wordsThymus_intRegul <- pubmed_parse˜_abstract(pubmedSearch_SSFintRegul$Thymus)
# # Search kidney signature
# wordsKidney_intRegul <- pubmed_parse_abstract(pubmedSearch_SSFintRegul$Kidney)
# # Search liver signature
# wordsLiver_intRegul <- pubmed_parse_abstract(pubmedSearch_SSFintRegul$Liver)
# # Search bone marrow signature
# wordsBoneMarrow_intRegul <- pubmed_parse_abstract(pubmedSearch_SSFintRegul$BoneMarrow)
# # Search Lung signature
# wordsLung_intRegul <- pubmed_parse_abstract(pubmedSearch_SSFintRegul$Lung)
# # Search SmallIntestine signature
# wordsSmallIntestine_intRegul <- pubmed_parse_abstract(pubmedSearch_SSFintRegul$SmallIntestine)
# # Search Spleen signature
# wordsSpleen_intRegul <- pubmed_parse_abstract(pubmedSearch_SSFintRegul$Spleen)
# # Search Testis signature
# wordsTestis_intRegul <- pubmed_parse_abstract(pubmedSearch_SSFintRegul$Testis)
# 
# words_df_list <- list(
#   Brain = wordsBrain_intRegul,
#   Thymus = wordsThymus_intRegul,
#   Kidney = wordsKidney_intRegul,
#   Liver = wordsLiver_intRegul,
#   BoneMarrow = wordsBoneMarrow_intRegul,
#   Lung = wordsLung_intRegul,
#   SmallIntestine = wordsSmallIntestine_intRegul,
#   Spleen = wordsSpleen_intRegul,
#   Testis = wordsTestis_intRegul
# )

# save it in html
library(htmlwidgets)
library(webshot)
# webshot::install_phantomjs()


lapply(names(words_SSFintRegul), function(id){
  words_df <- words_SSFintRegul[[id]]
  df <- words_df %>% dplyr::filter(!word %in% c("cell", "cells", "gene", "genes", "expression"))
  my_graph <- wordcloud2(data=df, shape="circle")
  saveWidget(my_graph, "tmp/tmp.html", selfcontained = FALSE)
  # and in png or pdf
  p <- file.path(dirname(params$p_nmfregulonAUC), paste0("wordCloudSSF_", id, ".pdf"))
  size = 1500
  webshot("tmp/tmp.html",p, delay =5, vwidth = size, vheight=size)
})

lapply(words_SSFintRegul, head)


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
