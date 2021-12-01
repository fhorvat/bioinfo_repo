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

library(ggrepel)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# get ggplot colors
gg_color_hue <- function(n){
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

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
  tibble(kmer_n = c(rep("6mer", 3), 
                    rep("7mer", 2)), 
         kmer_seq = c(c("CTGCTA", "GCTGCT", "TGCTGC"), 
                      c("GCTGCTA", "TGCTGCT")), 
         label = c(c("1-6", "2-7", "3-8"), 
                   c("1-7", "2-8")))

### plot sylamer output
purrr::map(names(syl_tb_list), function(kmer){
  
  # filter kmer table
  top_kmers <- 
    kmer_tb %>% 
    dplyr::filter(kmer_n == kmer) %$%
    kmer_seq
  
  # get syl out table
  kmer_tb_filt <- syl_tb_list[[kmer]]
  
  # get other kmers
  other_kmers <- 
    kmer_tb_filt %$%
    kmer %>% 
    unique(.) %>% 
    .[!(. %in% top_kmers)]
  
  # get labels
  kmer_tb_labels <-
    kmer_tb_filt %>% 
    dplyr::filter(kmer %in% top_kmers) %>%
    dplyr::group_by(kmer) %>% 
    dplyr::filter(abs(pvalue) == max(abs(pvalue))) %>% 
    dplyr::left_join(., kmer_tb, by = c("kmer" = "kmer_seq"))
  
  # set factor levels
  kmer_tb_levels <- 
    kmer_tb_filt %>% 
    dplyr::mutate(kmer = factor(kmer, levels = c(other_kmers, top_kmers)))
  
  # plot
  syl_plot <- 
    ggplot() +
    geom_line(data = kmer_tb_levels, aes(x = bin, y = pvalue, color = kmer, size = kmer)) +
    geom_label_repel(data = kmer_tb_labels, aes(label = str_c(label, ": ", kmer), x = bin, y = pvalue),
                     fontface = "bold", color = "black", box.padding = 0.35,
                     point.padding = 0.5, segment.color = "grey50") +
    scale_color_manual(values = c(rep("gray60", length(other_kmers)), 
                                  gg_color_hue(length(top_kmers)))) +
    scale_size_manual(values = c(rep(1, length(other_kmers)), 
                                 rep(2, length(top_kmers)))) +
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
    ggtitle(str_c("Sylamer enrichment - ", kmer, ", miR-15a seeds", sep = " ")) +
    ggsave(filename = file.path(outpath, sylamer_path[1] %>% basename %>% str_replace("6mer\\.txt$", str_c(kmer, ".seeds.png"))), 
           width = 12, height = 12)
  
  # return
  return(kmer)
  
})









