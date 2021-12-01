### INFO: 
### DATE: Wed Oct 24 14:39:04 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/other/abstract_wordclouds/whatsapp")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(tidytext)
library(wordcloud2)
library(htmlwidgets)
library(webshot)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# whatsapp chat path
whatsapp_path <- file.path(inpath, "WhatsApp_Chat_with_Ja_sam_tu.txt")

######################################################## READ DATA
# read whatsapp chat
whatsapp <- readr::read_lines(whatsapp_path)



######################################################## MAIN CODE
# tidy whatsapp messages 
whatsapp_df <- 
  tibble::tibble(raw = whatsapp) %>% 
  dplyr::mutate(date = str_extract(raw, "\\d+/\\d+/\\d+"), 
                time = str_extract(raw, "(?<=, )\\d{2}:\\d{2}"), 
                message = str_remove(raw, date) %>% 
                  str_remove(., time) %>% 
                  str_remove(., "^,  - "), 
                sender = str_extract(message, "^.*?(?=\\:)") %>%  
                  str_replace_all(., magrittr::set_names(c("Tea", "Iva", "Ivica", "Dubravka", "Filip", "Jasna", "Marko", "Zvonimir"), 
                                                         c("1 Love", "\\+385 91 595 1908", "\\+385 98 208 295", "Duda tele2", 
                                                           "Filip Horvat", "Jasna Grskovic", "Marko G", "Zvonac"))), 
                message = str_remove(message, "^.*?: ")) %>% 
  dplyr::filter(!is.na(sender)) %>% 
  dplyr::select(-raw)  
  
# # set pattern and replacment words
# pattern_words <- c("RNAs", "cells", "elavl2", "nanog", "oet", "pathways", "genes", "elements", "proteins", "embryos", "ltrs", "mescs", "sn", 
#                    "microRNA", "blastomeres", "lin28a", "transcriptional", "dcp1a", "lin28b", "nSN", "nsn", "retrotransposons", "tools", 
#                    "transcripts", "dcp2", "dna", "dicer", "genomes", "mice", "mammals", "mammalian", "developmental", "functional", 
#                    "oocytes")
# replacement_words <- c("RNA", "cell", "Elavl2", "Nanog", "OET", "pathway", "gene", "element", "protein", "embryo", "LTRs", "mESCs", "SN", 
#                        "miRNA", "blastomere", "Lin28A", "transcription", "Dcp1A", "Lin28B", "NSN", "NSN", "retrotransposon", "tool", 
#                        "transcript", "Dcp2", "DNA", "Dicer", "genome", "mouse", "mammal", "mammal", "development", "function", 
#                        "oocyte")

# split abstracts to words
whatsapp_tidy <- 
  tibble::tibble(message = whatsapp_df$message) %>% 
  tidytext::unnest_tokens(word, message, to_lower = T) %>% 
  # dplyr::anti_join(stop_words, by = "word") %>% 
  # dplyr::mutate(word = str_replace(word, "(?<!mate|alte|pate)rna", "RNA")) %>% 
  # dplyr::mutate(word = str_replace_all(word, magrittr::set_names(replacement_words, pattern_words))) %>% 
  # dplyr::filter(!str_detect(word, "^[:digit:]+$|^O$|review")) %>% 
  dplyr::count(word, sort = T) %>% 
  dplyr::slice(1:200) %>%
  dplyr::rename(freq = n) %>% 
  as.data.frame(.)

### plot text cloud
# create wordcloud
my_graph <- wordcloud2::wordcloud2(data = whatsapp_tidy, 
                                   size = 1.3, 
                                   # color = colorRampPalette(rev(c("lightblue", "cornflowerblue", "dodgerblue4", "deepskyblue4")))(20),
                                   color = colorRampPalette(rev(c("black", "orange", "red")))(20),
                                   ellipticity = 0.5, 
                                   rotateRatio = 0)

# save
htmlwidgets::saveWidget(my_graph, file.path(outpath, "whatsapp.jasamtu.html"), selfcontained = F)

# save as picture
webshot("whatsapp.jasamtu.html", "whatsapp.jasamtu.png", delay = 30, vwidth = 2000, vheight = 1000)

