rm(list = ls(all.names = TRUE))
# LOAD LIBRARIES
require(dplyr)
require(igraph)
require(LDAvis)
require(ldatuning)
require(NLP)
require(quanteda)
require(readtext)
require(rebus)
require(qdap)
require(SnowballC)
require(stringr)
require(slam)
require(svDialogs)
require(textclean)
require(textstem)
require(tidyverse)
require(tidytext)
require(topicmodels)
require(tm)
require(tsne)
#require(word2vec)

source("functions.R")


setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# SETTINGS
input_file <- "231-1_main"
output_dir <- paste("./output/", input_file, "/", sep="")
dir.create(output_dir)


# PREPROCESS INPUT
input <- readtext(paste("./input/", input_file, ".csv", sep=""), text_field = "text", docid_field = "ID")
lts_text <- input
lts_text$text <- replace_non_ascii(lts_text$text, remove.nonconverted = F)
lts_text <- as_tibble(lts_text)
lts_text$text <- lts_text$text %>%
  str_to_lower() %>%  # convert all the string to low alphabet
  str_replace_all(pattern = "\\bdont\\b ", replacement = "don't") %>%
  str_replace_all(pattern = "dk", replacement = " don't know") %>%
  replace_misspelling() %>% # replace misspellings
  replace_internet_slang() %>% # replace internet slang to normal words
  replace_contraction() %>% # replace contraction to their multi-word forms
  replace_word_elongation() %>% # replace informal writing with known semantic replacements
  str_to_lower() %>%  # convert all the string to low alphabet
  str_replace_all(pattern = "\\bnot ", replacement = "not") %>%
  str_replace_all(pattern = "\\bno ", replacement = "no") %>%
  removePunctuation() %>% # remove punctuation
  str_squish() %>% # reduces repeated whitespace inside a string.
  str_trim() # removes whitespace from start and end of string
custom_stop_words <- read.csv("custom_stop_words.csv", header = T)
#custom_stop_words <- as_tibble(custom_stop_words)
#full_stop_words <- bind_rows(stop_words, custom_stop_words)
full_stop_words <- as_tibble(custom_stop_words)
#tidy_lts <- lts_text %>% unnest_tokens(word,text, token="ngrams", n=2) %>% anti_join(full_stop_words)
tidy_lts <- lts_text %>% unnest_tokens(word,text) %>% anti_join(full_stop_words)
tidy_lts <- drop_na(tidy_lts, word)
tidy_lts <- mutate(tidy_lts, lemma = lemmatize_words(word))
tidy_lts <- mutate(tidy_lts, stem = stem_words(lemma))
keep_words <- tidy_lts %>% select(lemma)
tidy_review <- tidy_lts %>% select(doc_id,lemma) %>% inner_join(keep_words) %>% count(doc_id,lemma)


# STANDARDIZE SEMANTICALLY SIMILAR WORDS
words <- tidy_review[order(-tidy_review$n),] %>% select(lemma) %>% unique() %>% unlist()
tidy_lts_mod <- tidy_lts %>% select(doc_id, lemma)
checked_words <- c()
changed_words <- matrix(nrow=0, ncol=2)
for(word in words){
  syns <- synonyms(word, multiwords=F, return.list=F)
  found_syns <- intersect(checked_words, syns)
  if(length(found_syns) > 0){
    for(s in found_syns){
      yn <- askYesNo(paste0(word, " -> ", s))
      if(yn == T){
        tidy_lts_mod$lemma[tidy_lts_mod$lemma == word] <- s
        changed_words <- rbind(changed_words,c(word, s))
        break
      }
    }
    if(yn != T){
      checked_words <- append(checked_words, word)
    }
  }else{
    checked_words <- append(checked_words, word)
  }
}

write.csv(changed_words, file = paste(output_dir, input_file, "_changed_words.csv", sep=""), row.names = F)


changed_words <- read.csv(paste(output_dir, input_file, "_changed_words.csv", sep=""))
tidy_lts_mod <- tidy_lts %>% select(doc_id, lemma)
for(i in 1:nrow(changed_words)){
  tidy_lts_mod$lemma[tidy_lts_mod$lemma == changed_words[i,1]] <- changed_words[i,2]
}


# RETAIN RELEVANT WORDS TO REDUCE SPARSITY
tidy_lts_mod <- mutate(tidy_lts_mod, stem = stem_words(lemma))
tidy_lts_mod$stem[tidy_lts_mod$stem == "notsur"] <- "notknow"
tidy_lts_mod$stem[tidy_lts_mod$stem == "noidea"] <- "notknow"
tidy_lts_mod$stem[tidy_lts_mod$stem == "noclu"] <- "notknow"
tidy_lts_mod$stem[tidy_lts_mod$stem == "unsur"] <- "notknow"
tidy_review_mod <- tidy_lts_mod %>% select(doc_id, stem) %>% group_by(stem) %>% mutate(n = n())
tidy_review_mod <- tidy_review_mod %>% filter(n > 2) # eliminate words appearing less than 3 times
tidy_review_mod <- as_tibble(tidy_review_mod)


# FINAL REVIEW OF RELEVANT WORDS
dtm_review <- tidy_review_mod %>% cast_dtm(doc_id, stem, n) %>% as.matrix()
tdm_review <- t(dtm_review)
write.csv(tdm_review, file = paste(output_dir, input_file, "_tdm.csv", sep=""))


# GET REMAINING IDs
IDs <- colnames(tdm_review)


# SETTINGS
SEED <- 123456789
max_depth <- 2


# INITIALIZE
topics <- list()
documents <- as.data.frame(IDs)
colnames(documents) <- c("doc_id")


#RECURSIVE TOPIC CONSTRUCTION
rLDA <- function(tidy_review, SEED, max_depth, topics, documents, output_dir, cur_depth = 1, parent = "Topic ") {
  
  if(cur_depth > max_depth)
    return(list('topics'=topics, 'documents'=documents))
  
  # DTM for topic modeling
  dtm_review <- tidy_review %>% cast_dtm(doc_id, stem, n) %>% as.matrix()
  tdm_review <- t(dtm_review)
  IDs <- colnames(tdm_review)
  dtm_lts <- tidy_review %>% cast_dtm(doc_id, stem, n)
  print(dtm_lts)
  
  # Choose number of topics
  dir.create(paste0(output_dir,"plots/"))
  find_topic_num(dtm_lts, SEED, paste0(output_dir,"plots/",parent))
  num_topics <- as.numeric(dlgInput(paste0("Enter number of topics for ", parent))$res)
  if(num_topics < 2 || length(num_topics) == 0)
    return(list('topics'=topics, 'documents'=documents))
  
  # Topic modelling code
  myLDA <- LDA(dtm_lts, k = num_topics, control = list(seed = SEED, verbose = 1))
  #LDAvis::serVis(topicmodels2LDAvis(myLDA), out.dir = output_dir, open.browser = TRUE)
  
  # See topics
  topics_LDA <- terms(myLDA, threshold=0.1/cur_depth)
  names(topics_LDA) <- lapply(names(topics_LDA), function(x) paste0(parent, str_extract(x,"\\d+"), "."))
  print(topics_LDA)
  topics <- c(topics, topics_LDA)
  
  # Assign document topics
  documents_LDA <- as.data.frame(myLDA@gamma) 
  colnames(documents_LDA) <- names(topics_LDA)
  documents_LDA$doc_id <- IDs
  documents <- left_join(documents, documents_LDA)
  
  # Find subtopics
  for(i in names(topics_LDA)){
    if(cur_depth > 1) {
      documents[i] <- documents[i]*documents[parent]
    }
    docs <- documents_LDA$doc_id[documents_LDA[i] >= 1/num_topics]
    tidy_review_child <- tidy_review %>% filter(doc_id %in% docs) %>% filter(!stem %in% unlist(topics_LDA))
    tidy_review_child <- tidy_review_child %>% select(doc_id,stem) %>% group_by(stem) %>% mutate(n = n())
    ret <- rLDA(tidy_review_child, SEED, max_depth, topics, documents, output_dir, cur_depth + 1, i)
    topics <- ret$topics
    documents <- ret$documents
  }
  
  return(list('topics'=topics, 'documents'=documents))
  
}


t1 <- Sys.time()
ret <- rLDA(tidy_review_mod, SEED, max_depth, topics, documents, output_dir)
topics <- ret$topics
documents <- ret$documents


# OUTPUT TOPICS
indx <- sapply(topics, length)
topics_output <- as.data.frame(do.call(rbind,lapply(topics, `length<-`, max(indx))))
write.csv(topics_output, file = paste(output_dir, input_file, "_topics.csv", sep=""))
write.csv(documents, file = paste(output_dir, input_file, "_documents.csv", sep=""), row.names = F)


# FORMAT FINAL OUTPUT

# Number topics
num_base <- sum(grepl("^Topic \\d+\\.$", colnames(documents)))
num_leafs <- sum(grepl(strrep("\\d+\\.",max_depth), colnames(documents)))
num_total <- num_base + num_leafs

# Document topics
output_pt2 <- left_join(input, documents)
output_pt2 <- output_pt2 %>% mutate(topic1 = apply(.[4:(3+num_total)], 1, function(x) maxn(1,x,num_base)))
output_pt2 <- output_pt2 %>% mutate(topic2 = apply(.[4:(3+num_total)], 1, function(x) maxn(2,x,num_base)))
output_pt2 <- output_pt2 %>% select(doc_id, text, weight, topic1, topic2, everything())

# Topic weights
weights1 <- output_pt2 %>% group_by(topic1) %>% summarise(sum(weight)) %>% mutate(topic1 = replace(topic1, is.na(topic1), "NA"))
weights2 <- output_pt2 %>% group_by(topic2) %>% summarise(sum(weight)) %>% filter(!is.na(topic2))
weights_total <- merge(weights1, weights2, by.x="topic1", by.y="topic2", all=T) 
weights_total <- weights_total %>%
  mutate(weight = rowSums(weights_total[,c(2,3)], na.rm=T)) %>%
  mutate(percent = weight/sum(input$weight))
rownames(weights_total) <- weights_total[,1]
weights_total <- weights_total %>% select(weight, percent)
output_pt1 <- merge(weights_total, topics_output,by="row.names",all.x=TRUE)

# Write output
sink(paste(output_dir, input_file, "_output.csv", sep=""))
write.csv(output_pt1, row.names = F)
cat('\n')
write.csv(output_pt2, row.names = F)
sink()


t2 <- Sys.time()
t2-t1