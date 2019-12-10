### INFO: 
### DATE: Thu Aug 15 15:50:41 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/breeding")

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

library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# breeding table path
breeding_path <- file.path(inpath, "Breeding Performance Lnc1 190815.xlsx")

######################################################## READ DATA
# read breeding table
breeding_tb <- openxlsx::read.xlsx(breeding_path)

######################################################## MAIN CODE
# tidy breeding table
breeding_tidy <- 
  breeding_tb %>% 
  as_tibble(.) %>% 
  dplyr::select(male, female, litter_size = Total, mother_age = "Age.of.mother.in.days") %>% 
  dplyr::mutate(female_genotype = str_replace_all(female, c("\\+/\\+" = "WT", "\\+/-" = "Het", "-/-" = "Null")), 
                genotype = female_genotype) %>% 
  dplyr::mutate(mating_age = mother_age - 21) %>% 
  dplyr::filter(litter_size > 0) %>%  
  dplyr::mutate(genotype = factor(genotype, levels = c("WT", "Het", "Null"))) %>% 
  dplyr::select(genotype, mother_age = mating_age, litter_size)

# plot
breeding_plot <- 
  ggplot() +
  geom_point(data = breeding_tidy,
             mapping = aes(x = mother_age, y = litter_size, color = genotype), 
             size = 3.5) +
  scale_color_manual(values = c("#70ad47", "#bfbfbf", "black")) +
  scale_y_continuous(breaks = 0:max(breeding_tidy$litter_size)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none") +
  ggsave(filename = file.path(outpath, "litter_size.female_genotype.mating_age.scatter.png"), width = 10, height = 10)
