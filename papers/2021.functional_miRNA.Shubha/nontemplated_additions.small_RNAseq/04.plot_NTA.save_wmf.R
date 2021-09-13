### INFO: 
### DATE: Tue Aug 03 14:21:57 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/Svoboda/miRNA.Shubha/2021_paper.functional_oocyte_miRNA/maternal_miRNA_expression/nontemplated_additions/plot")
# setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/functional_oocyte_miRNA/maternal_miRNA_expression/nontemplated_additions/bam_subset")

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

# list of tables with values
softclip_tb_path <- list.files(inpath, "miRNA_3p_NTA\\.barplot.*\\.csv")
  
######################################################## READ DATA
# read tables
mirna_animals <- 
  purrr::map(softclip_tb_path, readr::read_csv) %>% 
  set_names(., str_extract(softclip_tb_path, "mouse.mm10.GarciaLopez_2015_RNA_GSE59254|cow.bosTau9|pig.susScr11"))

######################################################## MAIN CODE
# count sequences post 3' end for each animal
mirna_animals_counts <- purrr::map(names(mirna_animals), function(animal){
  
  # get the table into long format
  mirna_poly <- 
    mirna_animals[[animal]] %>% 
    tidyr::pivot_longer(cols = -c(gene_name, mirna_count), values_to = "nta_count", names_to = "softclip_seq") %>% 
    dplyr::mutate(nucleotide = str_extract(softclip_seq, "A|U|C|G") %>% factor(., levels = c("A", "U", "G", "C")), 
                  percentage = 100 * (nta_count / mirna_count), 
                  gene_name = factor(gene_name, levels = unique(gene_name)), 
                  softclip_seq = factor(softclip_seq, levels = rev(unique(softclip_seq)))) %>% 
    dplyr::select(gene_name, nucleotide, softclip_seq, percentage, nta_count, mirna_count)
  
  # plot
  nta_barplot <- 
    ggplot() + 
    geom_bar(data = mirna_poly, 
             mapping = aes(x = nucleotide, y = percentage, fill = softclip_seq), 
             width = 1, stat = "identity", position = "stack", color = "black") + 
    facet_wrap(~ gene_name, nrow = 1, strip.position = "bottom") +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(0, 60)) + 
    scale_fill_manual(values = c("A" = "#993333", "AA" = "#B76D6D", "AAA" = "#D4A8A8", "AAAA" = "#F3E3E3", 
                                 "U" = "#4166AF", "UU" = "#738EC5", "UUU" = "#A6B7DC", "UUUU" = "#D9E0F3", 
                                 "G" = "#7B6BB2", "GG" = "#9D91C6", "GGG" = "#BFB7DA", "GGGG" = "#E2DEEF",
                                 "C" = "#6F8940", "CC" = "#94A771", "CCC" = "#BAC5A3", "CCCC" = "#E0E4D5"), 
                      breaks = c("A", "AA", "AAA", "AAAA",
                                 "U", "UU", "UUU", "UUUU",
                                 "G", "GG", "GGG", "GGGG",
                                 "C", "CC", "CCC", "CCCC"),
                      labels = c("A", "2A", "3A", "4A",
                                 "U", "2U", "3U", "4U",
                                 "G", "2G", "3G", "4G",
                                 "C", "2C", "3C", "4C"),
                      drop = F) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(), 
          legend.title = element_blank()) +
    # theme(panel.border = element_blank()) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom", 
          strip.background = element_blank())
  
  # save
  ggsave(filename = file.path(outpath, str_c("miRNA_3p_NTA.barplot", animal, "wmf", sep = ".")),
         plot = nta_barplot,
         width = 12,
         height = 10,
         limitsize = FALSE, 
         device = "wmf")
  
  # return
  return(animal)
  
})


