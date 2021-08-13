#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other) in chunks
### DATE: Mon Mar 04 14:51:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/read_length_distribution/datasets/hamster_testis_Mov10l.8.5dpp.run_2.small_RNAseq")

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

# library size path
library_size_path <- file.path(inpath, "bam_files", "library_sizes.txt")

######################################################## READ DATA
# read class files
read_class_list <- 
  purrr::map(read_class_path, readr::read_csv) %>% 
  dplyr::bind_rows(.) 

# library sizes table
library_sizes <- readr::read_delim(library_size_path, delim = "\t",col_names = c("sample_id", "library_size")) 

######################################################## MAIN CODE
# get library sizes for 19-32nt reads
library_sizes %<>% 
  dplyr::filter(str_detect(sample_id, ".*\\.19to32nt$")) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt$"))

# set feature classes
class_hier <- 
  tibble(read_group = c("rRNA", "tRNA", 
                        "miRNA",
                        "SINE", "LINE", "LTR", "other_repeat",
                        "protein_coding", "other_ensembl", 
                        "not_annotated"), 
         class = c("rRNA", "tRNA",
                   "other", 
                   "retrotransposon", "retrotransposon", "retrotransposon", "repeat_other", 
                   "mRNA", "other", 
                   "other"))

# filter table
read_class_tb <- 
  read_class_list %>% 
  dplyr::filter(read_width >= 18, read_width <= 32) %>% 
  dplyr::left_join(., class_hier, by = "read_group") %>% 
  dplyr::group_by(sample_id, class, read_width) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., library_sizes, by = "sample_id") %>% 
  dplyr::mutate(rpm = (count / round(library_size / 1E6, 6))) %>% 
  dplyr::mutate(class = factor(class, levels = c("rRNA", "tRNA", "other", "mRNA", "repeat_other", "retrotransposon")))

# pivot to wide and save
read_class_sample <-
  read_class_tb %>%
  tidyr::pivot_wider(id_cols = c(sample_id, class), names_from = read_width, values_from = rpm) %>%
  dplyr::arrange(sample_id, class) %T>%
  readr::write_csv(., file.path(outpath, str_c("hamster_testis_Mov10l.8.5dpp.run_2.small_RNAseq",
                                               "read_class.RPM.20210706.csv", sep = ".")))

# get mean value per genotype, pivot to wide and save
read_class_mean <-
  read_class_tb %>%
  dplyr::mutate(genotype = str_extract(sample_id, "Mov10l1_KO|Mov10l1_WT")) %>%
  dplyr::group_by(genotype, class, read_width) %>%
  dplyr::summarise(rpm = mean(rpm)) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(genotype = factor(genotype, levels = c("Mov10l1_WT", "Mov10l1_KO"))) %>%
  tidyr::pivot_wider(id_cols = c(genotype, class), names_from = read_width, values_from = rpm) %>%
  dplyr::arrange(genotype, class) %T>%
  readr::write_csv(., file.path(outpath, str_c("hamster_testis_Mov10l.8.5dpp.run_2.small_RNAseq",
                                               "read_class.mean_RPM.20210706.csv", sep = ".")))


# get y-axis limit
y_lim <- 
  read_class_mean %>% 
  tidyr::pivot_longer(cols = -c(genotype, class), names_to = "read_width", values_to = "rpm") %>% 
  dplyr::group_by(genotype, read_width) %>% 
  dplyr::summarise(rpm = sum(rpm)) %$% 
  rpm %>% 
  max(.) %>% 
  ceiling(.)
  
# plot barplot separate for each genotype
purrr::map(unique(read_class_mean$genotype), function(genotype_name){
  
  # filter
  read_class_sample <- 
    read_class_mean %>% 
    dplyr::filter(genotype == genotype_name) %>% 
    dplyr::select(-genotype) %>% 
    tidyr::pivot_longer(cols = -class, names_to = "read_width", values_to = "rpm")
  
  # plot for each sample
  read_class_barplot <- 
    ggplot() + 
    geom_bar(data = read_class_sample, 
             mapping = aes(x = read_width, y = rpm, fill = class), width = 0.8, stat = "identity") + 
    scale_x_discrete(labels = as.character(c(18:32))) +
    scale_y_continuous(limits = c(0, 190)) +
    scale_fill_manual(values = c("rRNA" = "#F2F2F2", "tRNA" = "#BFBFBF", "other" = "#4472C4",
                                 "mRNA" = "#FFC000", "repeat_other" = "#ED7D31", "retrotransposon" = "#FF0000")) +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")
  
  # save
  ggsave(plot = read_class_barplot, filename = file.path(outpath, str_c(genotype_name, 
                                                                        "read_class.mean_RPM.20210706.barplot.png", 
                                                                        sep = ".")), 
         width = 12, height = 10)
  
  # return
  return(genotype_name)
  
})


