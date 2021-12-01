### INFO: counts hexamers in 3pUTRs 
### DATE: Thu Sep 16 01:31:05 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/functional_oocyte_miRNA/3pUTR_hexamers/hexamer_composition")

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

# hexamer count path
hexamer_path_list <- list.files(inpath, "\\.6mer_counts\\.csv", full.names = T)

######################################################## READ DATA
# read hexamer counts
hexamer_counts_tb <- purrr::map(hexamer_path_list, function(path){
  
  # read, add dataset
  readr::read_csv(path) %>% 
    dplyr::mutate(dataset = basename(path) %>% str_remove(., "\\.6mer_counts\\.csv"), 
                  animal = str_extract(dataset, "mouse|pig"), 
                  tissue = str_remove(dataset, "\\..*") %>% str_remove(., "^s_"))
  
}) %>%
  dplyr::bind_rows(.)
  
######################################################## MAIN CODE
# plot datasets
hexamer_plot_tb <- 
  hexamer_counts_tb %>% 
  tidyr::unite(col = "animal_tissue", animal, tissue) %>% 
  dplyr::group_by(animal_tissue) %>% 
  dplyr::mutate(hex_freq = count / sum(count))

# plot as points
kmer_points <- 
  ggplot(hexamer_plot_tb, aes(x = pos, y = hex_freq), fill = "black") + 
  facet_wrap(~ animal_tissue, ncol = 1, strip.position = "top") +
  geom_point() +
  ylab("") + 
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save 
ggsave(filename = file.path(outpath, str_c("all_tissue", "6mer_counts", "png", sep = ".")), 
       plot = kmer_points, 
       width = 10, height = 20)


### plot separately with same scale
# get y-scale
y_scale <- 
  hexamer_plot_tb %$% 
  hex_freq %>% 
  max(.)

# set manually
y_scale <- 0.0035
  
# split by dataset
hexamer_plot_list <- split(hexamer_plot_tb, hexamer_plot_tb$animal_tissue)

# plot
purrr::map(names(hexamer_plot_list), function(animal_tissue){
  
  # get the table
  hexamer_plot_tb <- 
    hexamer_plot_list[[animal_tissue]] %>% 
    dplyr::select(kmer, count, kmer_freq = hex_freq, pos, animal_tissue) %>% 
    dplyr::arrange(desc(kmer_freq))
  
  # # save
  # readr::write_csv(hexamer_plot_tb, str_c("freq_6mer", animal_tissue, "csv", sep = "."))
  
  # highlight some miRNA seeds
  mirna_seeds <- c("AATAAA", "GCACTT", "CACTCC", "CATTCC")

  # filter the table
  mirna_seeds_tb <- 
    hexamer_plot_tb %>% 
    dplyr::filter(kmer %in% mirna_seeds)
  
  # plot as points
  kmer_points <- 
    ggplot() + 
    geom_point(data = hexamer_plot_tb, mapping = aes(x = pos, y = kmer_freq), fill = "black") +
    geom_point(data = mirna_seeds_tb, mapping = aes(x = pos, y = kmer_freq), color = "red", size = 5) +
    # geom_text(data = mirna_seeds_tb, mapping = aes(x = pos, y = kmer_freq, label = kmer),
    #           check_overlap = TRUE, size = 3, hjust = 0, vjust = 1.5,
    #           colour = "black", fontface = "plain") + 
    geom_vline(xintercept = c(1, 1024, 2048, 3072, 4096), linetype = "dotted") + 
    scale_x_continuous(breaks = mirna_seeds_tb$pos,
                       labels = mirna_seeds_tb$kmer) +
    scale_y_continuous(limits = c(0, y_scale)) + 
    xlab("6mer rank sorted") +
    ylab("6mer frequency") + 
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 15, vjust = 0.5, angle = 90),
          axis.text.y = element_text(size = 15, vjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save 
  ggsave(filename = file.path(outpath, str_c("freq_6mer", animal_tissue, "png", sep = ".")), 
         plot = kmer_points, 
         width = 10, height = 10)
  
  # return
  return(animal_tissue)
  
})
