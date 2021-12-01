### INFO: summarize read classes
### DATE: Mon Mar 04 14:51:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/datasets/rodent_oocytes.small_RNAseq.2021_Sep/Analysis/rat/read_width_classification")

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

library(ggrepel)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# list read class files
read_class_path <- list.files(inpath, "\\.read_class\\.widths\\.csv", full.names = T, recursive = F)

######################################################## READ DATA
# read class files
read_class_list <- 
  purrr::map(read_class_path, readr::read_csv) %>% 
  dplyr::bind_rows(.) 

######################################################## MAIN CODE
# set feature classes
class_hier <- tibble(read_group = c("miRNA.mature.sense", "miRNA.other.sense",
                                    "protein_coding.sense",
                                    "rRNA", 
                                    "SINE", "LINE", "LTR", "other_repeat",
                                    "annotated_pseudogene", "other", "not_annotated"), 
                     class = c("miRNA", "miRNA", 
                               "mRNA", 
                               "other_mapped", 
                               "repeats", "repeats", "repeats", "repeats", 
                               "other_mapped", "other_mapped", "other_mapped"))

# get table
read_class_tb <- 
  read_class_list %>% 
  dplyr::filter(read_width >= 18, read_width <= 32) %>% 
  dplyr::left_join(., class_hier, by = "read_group")


### get rough classes
# prepare table for plot
read_class_rough <- 
  read_class_tb %>%
  dplyr::group_by(sample_id, class, read_width) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::group_by(sample_id, read_width) %>% 
  dplyr::mutate(total_count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(percentage = round((count / total_count), 3) * 100, 
                class = factor(class, levels = rev(unique(class_hier$class)))) %>% 
  dplyr::arrange(sample_id, read_width, desc(class)) %>% 
  dplyr::group_by(sample_id, read_width) %>% 
  dplyr::mutate(pos = cumsum(percentage) - (0.5 * percentage)) %>%
  dplyr::ungroup(.) %>% 
  dplyr::mutate(read_width = as.character(read_width) %>% factor(., levels = as.character(18:32)))

# plot for each sample
read_class_barplot <- 
  ggplot() + 
  geom_bar(data = read_class_rough, 
           mapping = aes(x = read_width, y = percentage, fill = class), 
           width = 0.8, stat = "identity") + 
  geom_label(data = read_class_rough,
             mapping = aes(x = read_width, y = pos, label = str_c(percentage, "%")),
             size = 4, stat = "identity", color = "black", fill = "white") +
  scale_fill_manual(values = c(miRNA = "#70ad47",
                               mRNA = "#ffc000",
                               repeats = "#ff0000",
                               other_mapped = "#000000")) +
  facet_grid(rows = vars(sample_id)) +
  scale_x_discrete(labels = as.character(c(18:32))) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right")

# save
ggsave(plot = read_class_barplot, 
       filename = file.path(outpath, "read_length.rough_classes.barplot.20211005.pdf"), 
       width = 12, height = 12)


### get precise classes
# prepare table for plot
read_class_precise <- 
  read_class_tb %>%
  dplyr::group_by(sample_id, read_group, read_width) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::group_by(sample_id, read_width) %>% 
  dplyr::mutate(total_count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(percentage = round((count / total_count), 3) * 100, 
                read_group = factor(read_group, levels = rev(unique(class_hier$read_group)))) %>% 
  dplyr::arrange(sample_id, read_width, desc(read_group)) %>% 
  dplyr::group_by(sample_id, read_width) %>% 
  dplyr::mutate(pos = cumsum(percentage) - (0.5 * percentage)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(read_width = as.character(read_width) %>% factor(., levels = as.character(18:32)))

# plot for each sample
read_class_barplot <- 
  ggplot() + 
  geom_bar(data = read_class_precise, 
           mapping = aes(x = read_width, y = percentage, fill = read_group), 
           width = 0.95, stat = "identity", color = "black") + 
  geom_label(data = read_class_precise %>% dplyr::filter(percentage > 2),
             mapping = aes(x = read_width, y = pos, label = str_c(percentage, "%")),
             size = 2, stat = "identity", color = "black", fill = "white") +
  # scale_fill_viridis_d(guide = guide_legend(reverse = TRUE)) +
  facet_grid(rows = vars(sample_id)) +
  scale_x_discrete(labels = as.character(c(18:32))) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right")

# save
ggsave(plot = read_class_barplot, 
       filename = file.path(outpath, "read_length.precise_classes.barplot.20211005.pdf"), 
       width = 12, height = 12)


