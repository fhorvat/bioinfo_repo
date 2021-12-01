#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other) in chunks
### DATE: Mon Mar 04 14:51:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/small_RNAseq_read_distribution/datasets/hamster_oocyte_Mov10l.smallRNAseq")

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
class_hier <- 
  tibble(read_group = c("miRNA",
                        "rRNA", "tRNA", 
                        "SINE", "LINE", "LTR", "other_repeat",
                        "protein_coding", "other_ensembl", 
                        "not_annotated"), 
         class = c("miRNA", 
                   "rRNA", "tRNA", 
                   "repeats", "repeats", "repeats", "repeats", 
                   "mRNA", "other_mapped", 
                   "other_mapped")) %>% 
  dplyr::mutate(class = read_group)

# filter table
read_class_tb <- 
  read_class_list %>% 
  dplyr::filter(str_detect(sample_id, "ovary")) %>% 
  dplyr::filter(!str_detect(sample_id, "19to32nt")) %>% 
  dplyr::left_join(., class_hier, by = "read_group") %>% 
  dplyr::group_by(sample_id, class, read_width) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(!(class %in% c("rRNA", "tRNA"))) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::mutate(total_count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(count_fract = (count / total_count)) %>% 
  dplyr::filter(read_width >= 18, read_width <= 32) %>% 
  dplyr::mutate(read_width = as.character(read_width) %>% factor(., levels = as.character(18:32)))

# plot for each sample
read_class_barplot <- 
  ggplot() + 
  geom_bar(data = read_class_tb, 
           mapping = aes(x = read_width, y = count_fract, fill = class), width = 0.8, stat = "identity") + 
  facet_grid(rows = vars(sample_id)) + 
  scale_x_discrete(labels = as.character(c(18:32))) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

# save
ggsave(plot = read_class_barplot, filename = file.path(outpath, "read_length.classes.18to32nt_reads.barplot.20210513.png"), width = 10, height = 20)

