### INFO: 
### DATE: Wed Oct 24 14:39:04 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/other/abstract_wordclouds")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(RefManageR)
library(tidytext)
library(wordcloud2)
library(htmlwidgets)
library(webshot)
library(RColorBrewer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# bibtex file path
bib_path <- file.path(inpath, "savedrecs.bib")

######################################################## READ DATA
# read bibtex
bib <- RefManageR::ReadBib(bib_path)

######################################################## MAIN CODE
# get abstracts and tidy them
abstracts <- 
  bib %>% 
  as.data.frame(.) %>% 
  as.tibble(.) %$% 
  abstract %>% 
  stringr::str_remove_all(., "\\n|\\{|\\}") %>% 
  stringr::str_remove_all(., "\\(C\\).*") %>% 
  stringr::str_squish(.) %>% 
  .[!is.na(.)]

# set pattern and replacment words
pattern_words <- c("RNAs", "cells", "elavl2", "nanog", "oet", "pathways", "genes", "elements", "proteins", "embryos", "ltrs", "mescs", "sn", 
                   "microRNA", "blastomeres", "lin28a", "transcriptional", "dcp1a", "lin28b", "nSN", "nsn", "retrotransposons", "tools", 
                   "transcripts", "dcp2", "dna", "dicer", "genomes", "mice", "mammals", "mammalian", "developmental", "functional", 
                   "oocytes")
replacement_words <- c("RNA", "cell", "Elavl2", "Nanog", "OET", "pathway", "gene", "element", "protein", "embryo", "LTRs", "mESCs", "SN", 
                       "miRNA", "blastomere", "Lin28A", "transcription", "Dcp1A", "Lin28B", "NSN", "NSN", "retrotransposon", "tool", 
                       "transcript", "Dcp2", "DNA", "Dicer", "genome", "mouse", "mammal", "mammal", "development", "function", 
                       "oocyte")

# split abstracts to words
abstract_tidy <- 
  tibble::tibble(abstract_text = abstracts) %>% 
  tidytext::unnest_tokens(word, abstract_text, to_lower = T) %>% 
  dplyr::anti_join(stop_words, by = "word") %>% 
  dplyr::mutate(word = str_replace(word, "(?<!mate|alte|pate)rna", "RNA")) %>% 
  dplyr::mutate(word = str_replace_all(word, magrittr::set_names(replacement_words, pattern_words))) %>% 
  dplyr::filter(!str_detect(word, "^[:digit:]+$|^O$|review")) %>% 
  dplyr::count(word, sort = T) %>% 
  dplyr::slice(1:100) %>%
  dplyr::rename(freq = n) %>% 
  as.data.frame(.)

### plot text cloud
# create wordcloud
my_graph <- wordcloud2::wordcloud2(data = abstract_tidy, 
                                   size = 1.3, 
                                   color = colorRampPalette(rev(c("lightblue", "cornflowerblue", "dodgerblue4", "deepskyblue4")))(20),
                                   # color = colorRampPalette(rev(c("black", "orange", "red")))(20),
                                   ellipticity = 0.5, 
                                   rotateRatio = 0)

# save
htmlwidgets::saveWidget(my_graph, file.path(outpath, "abstracts.2013to2018.top100words.blue.html"), selfcontained = F)

# save as picture
webshot("abstracts.2013to2018.top100words.blue.html", "abstracts.2013to2018.top100words.blue.png", delay = 30, vwidth = 2000, vheight = 1000)

