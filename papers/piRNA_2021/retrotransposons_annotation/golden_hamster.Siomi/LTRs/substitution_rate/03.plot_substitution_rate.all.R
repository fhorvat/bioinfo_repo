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

# table with LTRs info path
ltr_tb_path <- file.path(inpath, "LTRs.random_200_per_repName.200730.csv")

# chosen LTR classes path
classes_tb_path <- file.path(inpath, "..", "all_LTR_classes 200730.xlsx")

# substitution rate .RDS files path
sub_rate_path <- file.path(inpath, "sub_rate.RDS_files")
sub_rate_path <- list.files(sub_rate_path, pattern = "^sub_rate\\..*\\.RDS$", full.names = T)

######################################################## READ DATA
# read table with LTRs info
ltr_tb <- readr::read_csv(ltr_tb_path)

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
# prepare table for plot - add category and repFamily
sub_rate_tb_all %<>% 
  dplyr::left_join(., classes_tb, by = "repName") %>%
  dplyr::mutate(repFamily = factor(repFamily, levels = c("ERVK", "ERV1", "ERVL", "ERVL-MaLR", "Gypsy", "other")), 
                type = factor(type, levels = c("LTR", "INT"))) %>% 
  dplyr::arrange(repFamily, type) %>% 
  dplyr::mutate(repName = factor(repName, levels = unique(.$repName)))

### plot separately ERVK and all other
# split to ERVK and other
sub_rate_list <- split(sub_rate_tb_all, sub_rate_tb_all$repFamily == "ERVK")
names(sub_rate_list) <- c("all_other", "ERVK")

# loop through list
purrr::map(names(sub_rate_list), function(name){
  
  # get one table
  sub_rate_tb <- sub_rate_list[[name]]
  
  # data statistics
  sub_rate_stats <- 
    sub_rate_tb %>% 
    group_by(repName) %>% 
    summarize(median = median(sub_rate), 
              sd = sd(sub_rate), 
              length = length(sub_rate)) %>% 
    dplyr::mutate(SEM = sd / sqrt(length))
  
  # labels for plot
  labels_tb <- 
    sub_rate_tb %>% 
    dplyr::group_by(repName) %>% 
    dplyr::summarise(count = n()) %>% 
    dplyr::filter(count > 0)
  
  # plot
  sub_rate_boxplot <- 
    ggplot() +
    stat_boxplot(data = sub_rate_tb, aes(x = repName, y = sub_rate), geom = "errorbar", ) +
    geom_boxplot(data = sub_rate_tb, aes(x = repName, y = sub_rate), outlier.colour = NULL, outlier.shape = NA) +
    geom_jitter(data = sub_rate_tb, aes(x = repName, y = sub_rate), alpha = 0.4, 
                colour = "black", size = 1, width = 0.1, height = 0, show.legend = F) +
    scale_x_discrete(labels = str_c(labels_tb$repName, 
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
         filename = file.path(outpath, str_c("LTRs.random_200_per_repName.200730", name, "boxplot.png", sep = ".")), 
         width = 40, 
         height = 10)
  
  
  # return
  return(name)
  
})
