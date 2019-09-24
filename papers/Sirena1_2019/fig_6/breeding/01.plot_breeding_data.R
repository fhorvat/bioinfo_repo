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
  dplyr::mutate(male_genotype = str_replace_all(male, c("\\+/\\+" = "WT", "\\+/-" = "Het", "-/-" = "Null")), 
                female_genotype = str_replace_all(female, c("\\+/\\+" = "WT", "\\+/-" = "Het", "-/-" = "Null")), 
                genotype = str_c(male_genotype, " x ", female_genotype)) %>% 
  dplyr::filter(litter_size > 0,
                mother_age <= 336) %>% 
  dplyr::mutate(mother_age_category = cut(mother_age, c(0, 112, 224, max(mother_age))), 
                mother_age_category = factor(mother_age_category, labels = c("<16 weeks", "16 - 32 weeks", "32 - 48 weeks"))) %>% 
  dplyr::mutate(mating_age = mother_age - 21,
                mating_age = cut(mating_age, c(0, 112, 224, max(mating_age))),
                mating_age = factor(mating_age, labels = c("<16 weeks", "16 - 32 weeks", "32 - 48 weeks"))) %>%
  dplyr::filter(genotype %in% c("WT x WT", "Het x Het", "Null x Null")) %>% 
  dplyr::mutate(genotype = factor(genotype, levels = c("WT x WT", "Het x Het", "Null x Null"))) %>% 
  dplyr::select(genotype, mother_age = mother_age_category, litter_size)

# # significance
# breeding_test <- 
#   breeding_tidy %>% 
#   dplyr::filter(mother_age == "<20 weeks")
# 
# t.test(breeding_test %>% dplyr::filter(genotype == "WT x WT") %$% litter_size, 
#        breeding_test %>% dplyr::filter(genotype == "Null x Null") %$% litter_size, 
#        alternative = "greater")
# 
# wilcox.test(breeding_test %>% dplyr::filter(genotype == "WT x WT") %$% litter_size, 
#             breeding_test %>% dplyr::filter(genotype == "Null x Null") %$% litter_size, 
#             alternative = "greater")

# data statistics
breeding_statistics <- 
  breeding_tidy %>% 
  group_by(genotype, mother_age) %>% 
  summarize(median = median(litter_size), 
            sd = sd(litter_size), 
            length = length(litter_size)) %>% 
  dplyr::mutate(SEM = sd / sqrt(length))

# plot
breeding_plot <- 
  ggplot() +
  geom_bar(data = breeding_statistics, 
           mapping = aes(x = mother_age, y = median, fill = genotype),
           width = 0.8, stat = "identity", alpha = 0.3, 
           position = position_dodge(width = 0.8)) +
  geom_dotplot(data = breeding_tidy, 
               mapping = aes(x = mother_age, y = litter_size, fill = genotype),
               binaxis = "y", binwidth = 1, dotsize	= 0.25, stackdir = "center", colour = "black", show.legend = F, stackratio = 0.5,
               position = position_dodge(width = 0.8, preserve = "total")) +
  stat_summary(data = breeding_tidy, 
               mapping = aes(x = mother_age, y = litter_size, color = genotype),
               fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.6, 
               position = position_dodge(width = 0.8)) +
  geom_errorbar(data = breeding_statistics,
                mapping = aes(x = mother_age, ymin = median - SEM, ymax = median + SEM, color = genotype), 
                width = 0.4, size = 1, 
                position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("#70ad47", "#bfbfbf", "black")) +
  scale_color_manual(values = c("#70ad47", "#bfbfbf", "black")) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(breaks = 0:max(breeding_tidy$litter_size)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none") +
  ggsave(filename = file.path(outpath, "litter_size.both_genotype.mating_age.plot.png"), width = 10, height = 10)
