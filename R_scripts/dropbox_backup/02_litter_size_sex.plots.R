### INFO: boxplot + jitter of litter size in two different genotype mice
### DATE: 06. 10. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Zuzka")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

######################################################## PATH VARIABLES
outpath <- getwd()

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
litter_df <- read_csv("Book3_ZL_littersizes_LacZKO_sex.csv")

######################################################## MAIN CODE
# clean input data
litter_df_clean <- 
  litter_df %>% 
  dplyr::select(1:3) %>% 
  tidyr::gather(genotype, litter_size) %>% 
  dplyr::filter(!is.na(litter_size), 
                litter_size != 0) %>% 
  tidyr::separate(genotype, into = c("genotype", "sex"), sep = " ") %>% 
  dplyr::mutate(sex = str_replace(sex, "females", "female"), 
                sex = replace(sex, is.na(sex), "HET"))

###### plots
# loop through sex
for(sex_filter in (c("male", "female"))){
  
  sex_filter <- "female"
  
  # filter data
  litter_df_clean_sex <- 
    litter_df_clean %>% 
    dplyr::filter(sex %in% c(sex_filter, "HET"))
  
  # data statistics
  litter_statistics <- 
    litter_df_clean_sex %>% 
    group_by(genotype) %>% 
    summarize(median = median(litter_size), 
              sd = sd(litter_size), 
              length = length(litter_size)) %>% 
    dplyr::mutate(SEM = sd / sqrt(length))
  
  # median + SD without box
  set.seed(1000)
  ggplot() +
    geom_jitter(data = litter_df_clean_sex, aes(x = genotype, y = litter_size, fill = genotype), 
                shape = 21, colour = "black", size = 3, width = 0.3, height = 0, show.legend = F) +
    stat_summary(data = litter_df_clean_sex, aes(x = genotype, y = litter_size), 
                 fun.y = median, fun.ymin = median, fun.ymax = median, 
                 geom = "crossbar", width = 0.5) +
    geom_errorbar(data = litter_statistics, aes(x = genotype, ymin = median - sd, ymax = median + sd), width = 0.3) +
    scale_fill_manual(values = c("red", "blue")) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_y_continuous(limits = c(-2, 15), breaks = 0:15) +
    xlab("genotype") +
    ylab("litter size") +
    ggtitle(label = str_c(sex_filter, " KO")) + 
    ggsave(filename = str_c("litter_size.", str_replace_all(sex_filter, " ", "_"), ".median.sd.2.pdf"), width = 10, height = 10)
  
  ### median + errorbars 
  # # median + SEM without box
  # ggplot() +
  #   geom_jitter(data = litter_df_clean_sex, aes(x = genotype, y = litter_size, fill = genotype), size = 3, width = 0.1, height = 0, show.legend = F) +
  #   stat_summary(data = litter_df_clean_sex, aes(x = genotype, y = litter_size), 
  #                fun.y = median, fun.ymin = median, fun.ymax = median, 
  #                geom = "crossbar", width = 0.5) +
  #   geom_errorbar(data = litter_statistics, aes(x = genotype, ymin = median - SEM, ymax = median + SEM), width = 0.3) +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1), 
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank()) +
  #   scale_y_continuous(limits = c(0, 15), breaks = 0:15) +
  #   xlab("genotype") +
  #   ylab("litter size") +
  #   ggtitle(label = str_c(sex_filter, " KO")) + 
  #   ggsave(filename = str_c("litter_size.", str_replace_all(sex_filter, " ", "_"), ".median.SEM.pdf"), width = 10, height = 10)
  
  # # median only
  # ggplot() +
  #   geom_jitter(data = litter_df_clean_sex, aes(x = genotype, y = litter_size, fill = genotype), size = 3, width = 0.1, height = 0, show.legend = F) +
  #   stat_summary(data = litter_df_clean_sex, aes(x = genotype, y = litter_size), 
  #                fun.y = median, fun.ymin = median, fun.ymax = median, 
  #                geom = "crossbar", width = 0.5) +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1), 
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank()) +
  #   scale_y_continuous(limits = c(0, 15), breaks = 0:15) +
  #   xlab("genotype") +
  #   ylab("litter size") +
  #   ggtitle(label = str_c(sex_filter, " KO")) + 
  #   ggsave(filename = str_c("litter_size.", str_replace_all(sex_filter, " ", "_"), ".median.pdf"), width = 10, height = 10)

}
