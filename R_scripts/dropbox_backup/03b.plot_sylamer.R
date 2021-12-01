### INFO: 
### DATE: Mon May 04 23:08:45 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/2020_paper/Freimer_microarray/Sylamer_analysis")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(purrr)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# Sylamer path
sylamer_path <- list.files(inpath, pattern = ".*\\.sylout\\.[6,7]mer\\.txt", full.names = T)

######################################################## READ DATA
# read sylamer out
syl_tb_list <- purrr::map(sylamer_path, function(path){
  
  # read and clean table
  syl_tb_tidy <- 
    readr::read_delim(file = path, delim = "\t") %>% 
    dplyr::rename(kmer = upper) %>% 
    tidyr::pivot_longer(., -kmer, names_to = "bin", values_to = "pvalue") %>% 
    dplyr::mutate(bin = as.numeric(bin))
  
  # return
  return(syl_tb_tidy)
  
}) %>%
  set_names(., str_extract(sylamer_path, "6mer|7mer"))

######################################################## MAIN CODE
# create kmer table 
kmer_tb <- 
  tibble(kmer_n = c("6mer", "7mer"), 
         kmer_seq = c("GCTGCT", "TGCTGCT"))

### plot sylamer output
purrr::map(names(syl_tb_list), function(kmer){
  
  # filter kmer table
  kmer_seq_filt <- 
    kmer_tb %>% 
    dplyr::filter(kmer_n == kmer) %$%
    kmer_seq
  
  # get syl out table
  kmer_tb_filt <- syl_tb_list[[kmer]]
  
  # filter kmer to highlight
  kmer_tb_highlight <- 
    kmer_tb_filt %>% 
    dplyr::filter(kmer == kmer_seq_filt)

  # plot
  syl_plot <- 
    ggplot() +
    geom_line(data = kmer_tb_filt, aes(x = bin, y = pvalue, color = kmer)) +
    geom_line(data = kmer_tb_highlight, aes(x = bin, y = pvalue), color = "red", size = 2) +
    scale_color_manual(values = rep("gray60", length(unique(kmer_tb_filt$kmer)))) +
    scale_x_continuous(breaks = seq(0, 16000, 2000), 
                       labels = seq(0, 16000, 2000)) +
    xlab("3'UTRs sorted from downregulated to upregulated") +
    ylab("log10(enrichment p-value)") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
          axis.title.y = element_text(size = 15, vjust = - 0.2),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(legend.position = "none") +
    ggtitle(str_c("Sylamer enrichment - ", kmer, kmer_seq_filt, sep = " ")) +
    ggsave(filename = file.path(outpath, sylamer_path[1] %>% basename %>% str_replace("6mer\\.txt$", str_c(kmer, ".png"))), 
           width = 12, height = 12)
  
  # return
  return(kmer)
  
})









