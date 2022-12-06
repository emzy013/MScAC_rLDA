# FUNCTION FOR STANDARDIZING SIMILAR WORDS WITH WORD2VEC
## Requires a word2vec model read in as the bin file
## Currently not used in main code
std_word2vec <- function(tidy_review, bin_file) {
  
  words <- tidy_review[order(-tidy_review$n),] %>% select(lemma) %>% unique()
  word_model <- read.word2vec(bin_file, normalize=TRUE)
  vocab <- summary(word_model, type='vocabulary')
  emb <- as.matrix(word_model)
  
  if(TRUE){
    vocab <- vocab %>% str_to_lower() %>% str_remove_all(pattern = "_.*") # Gensim models have this pattern, need to be removed
    rownames(emb) <- vocab
  }
  
  gibberish <- words$lemma[!(words$lemma %in% vocab)]
  words <- words %>% filter(!(lemma %in% gibberish))
  words_emb <- emb[words$lemma,]
  similarity_matrix <- word2vec_similarity(words_emb, words_emb, type="cosine")
  
}


# PLOT OPTIMAL NUMBER OF TOPICS
find_topic_num <- function(dtm_lts, SEED, name) {
  result <- FindTopicsNumber(
    dtm_lts,
    topics = seq(from = 2, to = 10, by = 1),
    metrics = c("Griffiths2004","CaoJuan2009", "Arun2010", "Deveaud2014"),
    method = "Gibbs",
    control = list(seed = SEED),
    mc.cores = 2L,
    verbose = TRUE
  )

  jpeg(paste0(name,".jpg"), width = 700)
  FindTopicsNumber_plot(result)
  dev.off()
}


# FUNCTIONS FOR FINDING MOST LIKELY TOPIC
maxn <- function(n, x, num_topics) {
  if(sum(!is.na(x)) == 0) 
    return(NA)
  
  x_leafs <- x[(num_topics+1):length(x)]
  
  if(sum(!is.na(x_leafs)) == 0) 
    x <- x[1:num_topics]
  else if(n == 1 && sort(x_leafs, na.last=F)[length(x_leafs)] < (1/3))
    x <- x[1:num_topics]
  else
    x <- x_leafs
  
  if(sort(x, na.last=F)[length(x)-(n-1)] >= (1/3))
    return(names(x)[order(x, decreasing = TRUE)[n]])
  else
    return(NA)
}
max2 <- function(x) sort(x, na.last=F)[length(x)-1]


# CREATE FUNCTION THAT CONVERTS LDA OUTPUT FOR LDAvis
# https://gist.github.com/trinker/477d7ae65ff6ca73cace
# Transform Model Output for Use with the LDAvis Package
topicmodels2LDAvis <- function(x){
  svd_tsne <- function(x) tsne(svd(x)$u)
  if (ncol(topicmodels::posterior(x)[["topics"]]) < 3) stop("The model must contain > 2 topics")
  mat <- x@wordassignments
  LDAvis::createJSON(
    phi = topicmodels::posterior(x)[["terms"]],
    theta = topicmodels::posterior(x)[["topics"]],
    vocab = colnames(topicmodels::posterior(x)[["terms"]]),
    #mds.method = jsPCA,
    mds.method = svd_tsne,
    doc.length = slam::row_sums(mat, na.rm = TRUE),
    term.frequency = slam::col_sums(mat, na.rm = TRUE),
    reorder.topics = FALSE
  )
}


