#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other) in chunks
### DATE: Mon Mar 04 14:51:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/small_RNAseq_read_distribution")

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
class_hier <- tibble(read_group = c("miRNA",
                                    "rRNA", "SINE", "LINE", "LTR", "other_repeat",
                                    "protein_coding", "other_ensembl", 
                                    "not_annotated"), 
                     class = c("miRNA", 
                               "other_mapped", "repeats", "repeats", "repeats", "repeats", 
                               "mRNA", "other_mapped", 
                               "other_mapped"))

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
  dplyr::select(sample_id, class, read_width, count_fract) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "(?<=Mov10l_)WT|KO|HET"), 
                age = str_extract(sample_id, "13dpp|21dpp|8.5_dpp") %>% str_replace(., "_", ""), 
                reseq = str_detect(sample_id, ".reseq"), 
                genotype = ifelse(reseq, str_c(genotype, ".reseq"), genotype)) %>% 
  dplyr::group_by(age, genotype, class, read_width) %>% 
  dplyr::summarise(count_fract = mean(count_fract)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(genotype == "WT") %>% 
  dplyr::mutate(read_width = as.character(read_width) %>% factor(., levels = as.character(19:32)), 
                age = factor(age, levels = c("8.5dpp", "13dpp", "21dpp")), 
                class = factor(class, levels = rev(c("miRNA", "repeats", "mRNA", "other_mapped")))) %>% 
  dplyr::arrange(age)

### save for each age
# split to list
read_class_tb_list <- split(read_class_tb, read_class_tb$age)

# loop through list
purrr::map(names(read_class_tb_list), function(testis_age){
  
  # get one table, transform to wide
  read_class_tb_list[[testis_age]] %>% 
    tidyr::pivot_wider(id_cols = class, names_from = read_width, values_from = count_fract, names_prefix = "r.") %>% 
    readr::write_csv(., file.path(outpath, str_c("read_length.classes.fractions", testis_age, "csv", sep = ".")))
  
})

# plot for each age
read_class_barplot <- 
  ggplot() + 
  geom_bar(data = read_class_tb, 
           mapping = aes(x = read_width, y = count_fract, fill = class), width = 0.8, stat = "identity") + 
  scale_fill_manual(values = c(miRNA = "#70ad47", mRNA = "#ffc000", repeats = "#ff0000", other_mapped = "#000000")) + 
  facet_grid(rows = vars(age)) + 
  scale_x_discrete(labels = as.character(c(19:32))) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

# save
ggsave(plot = read_class_barplot, filename = file.path(outpath, "read_length.classes.barplot.png"), width = 10, height = 10)



# # filter table
# read_class_tb <- 
#   read_class_list %>% 
#   dplyr::group_by(sample_id) %>% 
#   dplyr::mutate(total_count = sum(count)) %>% 
#   dplyr::ungroup(.) %>% 
#   dplyr::mutate(count_fract = (count / total_count)) %>% 
#   dplyr::select(sample_id, read_group, read_width, count_fract) %>% 
#   dplyr::mutate(genotype = str_extract(sample_id, "(?<=Mov10l_)WT|KO|HET"), 
#                 age = str_extract(sample_id, "13dpp|21dpp|8.5_dpp") %>% str_replace(., "_", ""), 
#                 reseq = str_detect(sample_id, ".reseq"), 
#                 genotype = ifelse(reseq, str_c(genotype, ".reseq"), genotype)) %>% 
#   dplyr::group_by(age, genotype, read_group, read_width) %>% 
#   dplyr::summarise(count_fract = mean(count_fract)) %>% 
#   dplyr::ungroup(.) %>% 
#   dplyr::filter(genotype == "WT") %>% 
#   dplyr::mutate(read_width = factor(read_width, levels = as.character(19:32)), 
#                 age = factor(age, levels = c("8.5dpp", "13dpp", "21dpp")), 
#                 read_group = factor(read_group, levels = rev(c("miRNA", "rRNA",
#                                                                "SINE", "LINE", "LTR", "other_repeat", 
#                                                                "protein_coding", "other_ensembl", 
#                                                                "not_annotated")))) %>% 
#   dplyr::arrange(age)
# 
# # plot for each age
# read_class_barplot <- 
#   ggplot() + 
#   geom_bar(data = read_class_tb, 
#            mapping = aes(x = read_width, y = count_fract, fill = read_group), width = 0.8, stat = "identity") + 
#   scale_fill_viridis_d() + 
#   facet_grid(rows = vars(age)) + 
#   scale_x_discrete(labels = as.character(c(19:32))) +
#   guides(fill = guide_legend(reverse = TRUE)) +
#   theme_bw() +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank()) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "bottom")
# 
# # save
# ggsave(plot = read_class_barplot, filename = file.path(outpath, "read_length.broad_classes.barplot.png"), width = 10, height = 10)
# 
# 
