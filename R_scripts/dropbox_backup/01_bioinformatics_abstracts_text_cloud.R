### INFO: get abstracts from NCBI and output top 100 words in wordcloud
### DATE: 24. 10. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()
# options(bitmapType = 'cairo')

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/NCBI_abstract_mine")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

library(RISmed)
library(tidytext)
library(wordcloud2)
library(htmlwidgets)

######################################################## PATH VARIABLES
outpath <- getwd()

######################################################## SOURCE FILES
# lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
# source(file.path(lib_path, "headt.R"))
# source(file.path(lib_path, "asdf.R"))

######################################################## FUNCTIONS

######################################################## READ DATA

######################################################## MAIN CODE
# set list of journals
journal_list <- c("Bioinformatics", "Cell", "Nature", "Science")

# get through the list and output top 100 words in all abstracts in 2016
for(journal in journal_list){
  
  # search pubmed for papers in Bioinformatics journal published in 2016, get abstracts
  query <- str_c("(\"", journal, "\"[Journal]) AND (\"2016/01/01\"[Date - Publication] : \"2016/12/31\"[Date - Publication])")
  cat(query, "\n")
  abstract_text <- 
    EUtilsSummary(query, type = "esearch", db = "pubmed", retmax = 30000) %>% 
    EUtilsGet(.) %>% 
    RISmed::AbstractText(.)  
  
  # set words to filter
  word_filter <- c("AVAILABILITY", "IMPLEMENTATION", "SUPPLEMENTARY", "INFORMATION", "MOTIVATION", 
                   "supplementary", "http", "CONTACT", "https", "github.com", "Supplementary")
  
  # tidy text
  abstract_tidy <- 
    unlist(abstract_text) %>% 
    .[stringr::str_length(.) > 0] %>% 
    tibble::tibble(abstract_text = .) %>% 
    tidytext::unnest_tokens(word, abstract_text, to_lower = FALSE) %>% 
    dplyr::filter(!(word %in% word_filter)) %>%
    dplyr::mutate(word = tolower(word)) %>% 
    anti_join(stop_words, by = "word") %>% 
    dplyr::mutate(word = replace(word, word == "methods", "method"), 
                  word = replace(word, word == "genes", "gene"), 
                  word = replace(word, word == "proteins", "protein"), 
                  word = replace(word, word == "networks", "network")) %>% 
    count(word, sort = T) %>% 
    dplyr::slice(1:100) %>% 
    dplyr::rename(freq = n) %>% 
    as.data.frame(.)
  
  # plot text cloud
  my_graph <- wordcloud2::wordcloud2(data = abstract_tidy, 
                                     size = 1.3, 
                                     color = colorRampPalette(rev(c("black", "red", "orange")))(20),
                                     ellipticity = 0.5, 
                                     rotateRatio = 0)
  saveWidget(my_graph, str_c(journal, "_top100_2016.html"), selfcontained = F)
  # webshot("tmp.html", "fig_1.pdf", delay = 5, vwidth = 1100, vheight = 600)
  
}

