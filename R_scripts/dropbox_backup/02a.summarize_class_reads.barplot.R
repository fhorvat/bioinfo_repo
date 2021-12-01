#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other) in chunks
### DATE: Mon Mar 04 14:51:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/antiviral_RNAi.Marcos/datasets/mouse_brain.TBEV.small_RNAseq.2021_Feb/Analysis/read_classification")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(data.table)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# list read class files
read_class_path <- list.files(inpath, "\\.read_class\\.widths\\.csv", full.names = T, recursive = T)

######################################################## READ DATA
# read class files
read_class_list <- 
  purrr::map(read_class_path, readr::read_csv) %>% 
  dplyr::bind_rows(.) 

######################################################## MAIN CODE
# set feature classes
class_hier <- tibble(read_group = c("miRNA.mature",
                                    "rRNA", "tRNA",
                                    "snoRNA", "snRNA", "protein_coding", 
                                    "not_annotated"), 
                     class = c("miRNA.mature",
                               "rRNA", "tRNA",
                               "snoRNA", "snRNA", "protein_coding", 
                               "not_annotated"))

# filter table
read_class_tb <- 
  read_class_list %>% 
  dplyr::left_join(., class_hier, by = "read_group") %>% 
  dplyr::group_by(sample_id, class, read_width) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::mutate(total_count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(count_fract = (count / total_count)) %>% 
  # dplyr::mutate(genotype = str_extract(sample_id, "WT|HET")) %>% 
  dplyr::mutate(read_width = as.character(read_width) %>% factor(., levels = as.character(8:80)), 
                class = factor(class, levels = rev(c("miRNA.mature",
                                                     "rRNA", "tRNA",
                                                     "snoRNA", "snRNA", "protein_coding", 
                                                     "not_annotated"))))

# plot as barplot
read_class_barplot <- 
  ggplot() + 
  geom_bar(data = read_class_tb, 
           mapping = aes(x = read_width, y = count_fract, fill = class), width = 0.8, stat = "identity") + 
  # scale_fill_manual(values = c(miRNA = "#70ad47", mRNA = "#ffc000", repeats = "#ff0000", other_mapped = "#000000")) + 
  facet_grid(rows = vars(sample_id)) + 
  scale_x_discrete(labels = as.character(c(8:80))) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

# save
ggsave(plot = read_class_barplot, filename = file.path(outpath, "read_length.classes.barplot.20210225.png"), width = 10, height = 10)
