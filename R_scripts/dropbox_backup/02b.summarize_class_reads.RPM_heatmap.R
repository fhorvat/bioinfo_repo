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

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/antiviral_RNAi.Marcos/datasets/mouse_brain.TBEV.small_RNAseq.2021_Feb/Data/Mapped/STAR_mm10"

# library size path
library_size_path <- file.path(mapped_path, "4_library_size", "library_sizes.txt")

######################################################## READ DATA
# read class files
read_class_list <- 
  purrr::map(read_class_path, readr::read_csv) %>% 
  dplyr::bind_rows(.) 

# read library size df
library_size_df <-  readr::read_delim(file = library_size_path, delim = "\t", col_names = c("sample_id", "library_size"))

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
  dplyr::left_join(., library_size_df, by = "sample_id") %>%
  dplyr::mutate(library_size = round((library_size / 1e6), 4),
                rpm = count / library_size) %>% 
  dplyr::filter(read_width >= 16, read_width <= 78) %>% 
  dplyr::mutate(read_width = as.character(read_width) %>% factor(., levels = as.character(16:78)), 
                class = factor(class, levels = rev(c("miRNA.mature",
                                                     "rRNA", "tRNA",
                                                     "snoRNA", "snRNA", "protein_coding", 
                                                     "not_annotated")))) %>% 
  dplyr::mutate(rpm = log10(rpm + 0.001), 
                count = log10(count + 1))

### plot for all samples
purrr::map(unique(read_class_tb$sample_id), function(sample_name){
  
  # create plot table
  plot_tb <- 
    read_class_tb %>% 
    dplyr::filter(sample_id == sample_name)
  
  # plot as heatmap ggplot2
  heat_plot <- 
    ggplot(plot_tb, aes(x = read_width, y = class)) + 
    geom_tile(aes(fill = count), colour = "grey45") + 
    # coord_equal() + 
    scale_fill_gradient2(low = "gray95", mid = "cyan", high = "darkgreen") +
    theme(axis.text.x = element_text(size = 12, angle = 90, face = "bold", colour = "grey25", vjust = 0.5, hjust = 0), 
          axis.text.y = element_text(size = 12, face = "bold", colour = "grey25"), 
          legend.title = element_text(size = 10, face = "bold"), 
          # legend.position = "rigth", 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA, colour = NA), 
          axis.ticks = element_blank()) + 
    labs(x = "", 
         y = "", 
         fill = "") + 
    ggtitle(sample_name) +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(limits = levels(plot_tb$class), position = "right")
  
  # save
  ggsave(plot = heat_plot, filename = file.path(outpath, str_c("read_length.classes.count_heatmap.", sample_name, ".20210225.png")), width = 15, height = 6)
         
  # return
  return(sample_name)
  
})

