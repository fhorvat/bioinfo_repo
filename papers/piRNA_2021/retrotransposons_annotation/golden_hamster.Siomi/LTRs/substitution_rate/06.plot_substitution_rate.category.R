### INFO: 
### DATE: Wed Jul 15 10:34:22 2020
### AUTHOR: fhorvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LTRs/substitution_rate/random_sampled_LTRs.all_classes")

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

# chosen LTR classes path
classes_tb_path <- file.path(inpath, "..", "all_LTR_classes 200730.xlsx")

# substitution rate .RDS files path
sub_rate_path <- file.path(inpath, "sub_rate.RDS_files")
sub_rate_path <- list.files(sub_rate_path, pattern = "^sub_rate\\..*\\.RDS$", full.names = T)

######################################################## READ DATA
# read table with chosen LTR classes
classes_tb <- openxlsx::read.xlsx(classes_tb_path)

# read substitution rate files, bind to one table
sub_rate_tb_all <- purrr::map(sub_rate_path, function(path){
  
  # get name
  fasta_name <- path %>% basename(.) %>% str_remove_all(., "^sub_rate\\.|\\.RDS$")
  
  # read from .RDS
  sub_rate_raw <- readRDS(path)
  
  # add to tibble
  sub_rate_raw_tb <- tibble(sub_rate = sub_rate_raw, 
                            repName = fasta_name)
  
}) %>% 
  bind_rows(.)

######################################################## MAIN CODE
### clean classes table
# select relevant columns
classes_tb %<>% 
  as_tibble(.) %>% 
  dplyr::select(repName, repFamily, category_I, type) %>% 
  dplyr::mutate(repFamily = replace(repFamily, is.na(repFamily), "other"))

### plot
# set subset of categories
category_sub <- c("LTR7x ERV1", "LTRISx", "MERx ERV1", "MLTRx ERV1", 
                  "MMERGLN", "MuLV", "RLTR1x ERV1", "RodERV21", 
                  "ERVBx", "ETn", "IAP", "MLTRx ERVK", 
                  "MMERVKx", "MYSERVx", "MuERV4", "RLTRx", 
                  "RLTR1x ERVK", "RLTR2x ERVK", "RLTR3x ERVK", "RLTR4x ERVK", 
                  "RMERx ERVK", "ERV3", "ERVL-x", "LTR6x ERVL", "LTR7x ERVL", 
                  "MERx ERVL", "MERVL", "RMERx ERVL", "MLT", 
                  "MT", "MT2", "ORR1", "Gypsy")

# prepare table for plot - add category, filter, subset
sub_rate_tb <- 
  sub_rate_tb_all %>% 
  dplyr::left_join(., classes_tb, by = "repName") %>% 
  dplyr::filter(category_I %in% category_sub) %>% 
  dplyr::group_by(category_I, type) %>% 
  dplyr::slice_sample(n = 200) %>% 
  dplyr::mutate(category_I = factor(category_I, levels = category_sub))

# data statistics
sub_rate_stats <- 
  sub_rate_tb %>% 
  group_by(category_I) %>% 
  summarize(median = median(sub_rate), 
            sd = sd(sub_rate), 
            length = length(sub_rate)) %>% 
  dplyr::mutate(SEM = sd / sqrt(length))

# labels for plot
labels_tb <- 
  sub_rate_tb %>% 
  dplyr::group_by(category_I) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::filter(count > 0)

# plot
sub_rate_boxplot <- 
  ggplot() +
  stat_boxplot(data = sub_rate_tb, aes(x = category_I, y = sub_rate), geom = "errorbar", size = 1.5) +
  geom_boxplot(data = sub_rate_tb, aes(x = category_I, y = sub_rate), outlier.colour = NULL, outlier.shape = NA, size = 1.5) +
  geom_jitter(data = sub_rate_tb, aes(x = category_I, y = sub_rate), alpha = 0.4, 
              colour = "black", size = 1.5, width = 0.1, height = 0, show.legend = F) +
  scale_x_discrete(labels = str_c(labels_tb$category_I, 
                                  " (",
                                  labels_tb$count, 
                                  ")"), 
                   drop = TRUE) + 
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank())


# save plot
ggsave(plot = sub_rate_boxplot, 
       filename = file.path(outpath, str_c("LTRs.random_200_per_category.200805", "boxplot.png", sep = ".")), 
       width = 20, 
       height = 10)
