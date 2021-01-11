### INFO: 
### DATE: Wed Jul 15 10:34:22 2020
### AUTHOR: fhorvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY

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



######################################################## READ DATA
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LINE1/substitution_rate/random_sampled_LINE1s")

# set common outpath
common_outpath <- getwd()

# classes table path
# classes_tb_path <- file.path(common_outpath, "LINE1s.random_200_per_repSubfamily.classes.merged.200807.csv")
classes_tb_path <- file.path(common_outpath, "LINE1s.random_200_per_repSubfamily.classes.200807.csv")

# read classes table
classes_tb <- 
  readr::read_csv(file = classes_tb_path) %>% 
  dplyr::mutate(repSubfamily = factor(repSubfamily, levels = unique(.$repSubfamily)))

######################################################## MAIN CODE
### list of animals to plot
# create list
animals_list <- c("golden_hamster.Siomi_assembly", "mouse.mm10", "rat.rn6", "chinese_hamster.CriGri_PICR")

### loop through animals
sub_rate_tb <- purrr::map(animals_list, function(animal){
  
  # set inpath 
  inpath <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation", 
                      animal, 
                      "LINE1/substitution_rate/random_sampled_LINE1s")

  
  # substitution rate .RDS files path
  sub_rate_path <- file.path(inpath, "sub_rate.RDS_files")
  sub_rate_path <- list.files(sub_rate_path, pattern = "^sub_rate\\..*\\.RDS$", full.names = T)
  
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
  
  # prepare table for plot 
  set.seed(1234)
  sub_rate_tb <- 
    sub_rate_tb_all %>% 
    dplyr::inner_join(., classes_tb, by = "repName") %>%
    dplyr::group_by(repSubfamily) %>% 
    dplyr::slice_sample(n = 200) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::mutate(repSubfamily = factor(repSubfamily, levels = unique(classes_tb$repSubfamily))) %>%
    dplyr::arrange(repSubfamily) %>% 
    dplyr::mutate(animal = animal)
  
  # return
  return(sub_rate_tb)
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(animal = str_remove(animal, "\\..*"), 
                animal = factor(animal, levels = c("golden_hamster", "chinese_hamster", "mouse", "rat")))

# data statistics
sub_rate_stats <- 
  sub_rate_tb %>% 
  group_by(animal, repSubfamily) %>% 
  summarize(median = median(sub_rate), 
            sd = sd(sub_rate), 
            length = length(sub_rate)) %>% 
  dplyr::mutate(SEM = sd / sqrt(length))

# labels for plot
labels_tb <- 
  sub_rate_tb %>% 
  dplyr::group_by(animal, repSubfamily) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::filter(count > 0) %>%
  dplyr::right_join(classes_tb %>% dplyr::select(repSubfamily) %>% unique(.), by = "repSubfamily") %>% 
  dplyr::mutate(count = replace(count, is.na(count), 0)) %>% 
  dplyr::arrange(repSubfamily)

# plot
sub_rate_boxplot <- 
  ggplot() +
  stat_boxplot(data = sub_rate_tb, aes(x = repSubfamily, y = sub_rate), geom = "errorbar", size = 1.5) +
  geom_boxplot(data = sub_rate_tb, aes(x = repSubfamily, y = sub_rate), outlier.colour = NULL, outlier.shape = NA, size = 1.5) +
  geom_jitter(data = sub_rate_tb, aes(x = repSubfamily, y = sub_rate), alpha = 0.4, 
              colour = "black", size = 1.5, width = 0.1, height = 0, show.legend = F) +
  facet_grid(rows = vars(animal)) +
  scale_x_discrete(drop = F) + 
  scale_y_continuous(limits = c(0, 0.5)) +
  theme_bw(base_size = 30) +
  theme(axis.text.x = element_text(size = 30, angle = 90, hjust = 1, vjust = 0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank())


# save plot to common path too
ggsave(plot = sub_rate_boxplot, 
       filename = file.path(common_outpath, str_c("LINE1s.random_200_per_repSubfamily.200805", "all_facet", "complete", "boxplot.png", sep = ".")), 
       width = 15, 
       height = 20)

