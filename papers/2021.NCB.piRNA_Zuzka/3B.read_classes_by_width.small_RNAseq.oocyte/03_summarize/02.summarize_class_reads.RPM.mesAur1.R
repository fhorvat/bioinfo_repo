#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other) in chunks
### DATE: Mon Mar 04 14:51:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/small_RNAseq_read_distribution/datasets/hamster_oocyte_Mov10l.deduplicated.smallRNAseq/mesAur1")

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
  dplyr::filter(read_width >= 18, read_width <= 32) %>% 
  dplyr::left_join(., class_hier, by = "read_group") %>% 
  dplyr::group_by(sample_id, class, read_width) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., library_sizes, by = "sample_id") %>% 
  dplyr::mutate(rpm = (count / round(library_size / 1E6, 6)),
                class = factor(class, levels = class_hier$class))

# save table separate for each sample
purrr::map(unique(read_class_tb$sample_id), function(sample_name){

  # filter, pivot to wide and save
  read_class_sample <-
    read_class_tb %>%
    dplyr::filter(sample_id == sample_name) %>%
    tidyr::pivot_wider(id_cols = c(sample_id, class), names_from = read_width, values_from = rpm) %>%
    dplyr::arrange(class) %T>%
    readr::write_csv(., file.path(outpath, str_c(sample_name, "read_class.deduplicated.RPM.20210517.csv", sep = ".")))

  # return
  return(sample_name)

})

# plot barplot separate for each sample
purrr::map(unique(read_class_tb$sample_id), function(sample_name){
  
  # filter
  read_class_sample <- 
    read_class_tb %>% 
    dplyr::filter(sample_id == sample_name)  
  
  # plot for each sample
  read_class_barplot <- 
    ggplot() + 
    geom_bar(data = read_class_tb, 
             mapping = aes(x = read_width, y = rpm, fill = class), width = 0.8, stat = "identity") + 
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
  ggsave(plot = read_class_barplot, filename = file.path(outpath, str_c(sample_name, "read_class.deduplicated.RPM.20210517.barplot.png")), 
         width = 10, height = 20)
  
  # return
  return(sample_name)
  
})
  

