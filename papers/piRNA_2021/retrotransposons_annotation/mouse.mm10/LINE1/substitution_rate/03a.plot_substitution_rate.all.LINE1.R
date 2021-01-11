### INFO: 
### DATE: Wed Jul 15 10:34:22 2020
### AUTHOR: fhorvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/mouse.mm10/LINE1/substitution_rate/random_sampled_LINE1s")

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

# table with LINE1s info path
line1_tb_path <- file.path(inpath, "LINE1s.random_200_per_repName.200806.csv")

# substitution rate .RDS files path
sub_rate_path <- file.path(inpath, "sub_rate.RDS_files")
sub_rate_path <- list.files(sub_rate_path, pattern = "^sub_rate\\..*\\.RDS$", full.names = T)

######################################################## READ DATA
# read table with LINE1s info
line1_tb <- readr::read_csv(line1_tb_path)

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
# ### list of retrotransposons to plot
# # create tibble
# classes_tb <- tibble(repName = c("IAPLTR3", "IAPLTR4", "IAP1-MM_LTR", "IAPLTR1a_Mm", "IAPLTR1_Mm", 
#                                  "IAPLTR2a", "IAPLTR2a2_Mm", "IAPLTR2b", "IAPLTR2_Mm", "IAPEY_LTR", 
#                                  "IAPEY2_LTR", "IAPEY3C_LTR", "IAPEY3_LTR", "IAPEY4_LTR", "IAPEY5_LTR", 
#                                  "IAPLTR3-int", "IAPLTR4_I", "IAP1-MM_I", "IAP-d-int", "IAPEY3-int", 
#                                  "IAPEY4_I", "IAPEY5_I", "IAPEy-int", "IAPEz-int", "IAPA_MM-int"))

### plot
# prepare table for plot 
sub_rate_tb <- 
  sub_rate_tb_all %>% 
  # dplyr::right_join(., classes_tb, by = "repName") %>%
  # dplyr::mutate(repName = factor(repName, levels = unique(classes_tb$repName))) %>% 
  dplyr::arrange(repName)

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
  stat_boxplot(data = sub_rate_tb, aes(x = repName, y = sub_rate), geom = "errorbar", size = 1.5) +
  geom_boxplot(data = sub_rate_tb, aes(x = repName, y = sub_rate), outlier.colour = NULL, outlier.shape = NA, size = 1.5) +
  geom_jitter(data = sub_rate_tb, aes(x = repName, y = sub_rate), alpha = 0.4, 
              colour = "black", size = 1.5, width = 0.1, height = 0, show.legend = F) +
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
       filename = file.path(outpath, str_c("LINE1s.random_200_per_repName.20200827", "boxplot.png", sep = ".")), 
       width = 40, 
       height = 10)

