### INFO: 
### DATE: Fri Aug 16 13:47:24 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/mitochondria_microscopy")

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

# mitochondria table path
mito_path <- file.path(inpath, "20190816 Lnc1 mitochondria cluster count tables.xlsx")

######################################################## READ DATA
# read mitochondria table
mito_tb_list <- 
  purrr::map(c("Mean", "Median"), function(sheet) openxlsx::read.xlsx(xlsxFile = mito_path, sheet = sheet) %>% as_tibble(.)) %>% 
  set_names(c("mean", "median"))

######################################################## MAIN CODE
# tidy tables, plot
purrr::map(c("mean", "median"), function(mito_stat){
  
  # gather + separate table
  mito_tb_tidy <- 
    mito_tb_list[[mito_stat]]  %>% 
    select(-1) %>% 
    tidyr::gather(key = "stage", value = "value") %>% 
    tidyr::separate(stage, into = c("stage", "genotype"), sep = "\\.") %>% 
    dplyr::filter(!is.na(value)) %>% 
    dplyr::mutate(genotype = factor(genotype, levels = c("WT", "Null")))
  
  # wilcox.test(mito_tb_tidy %>% dplyr::filter(stage == "GV") %>% dplyr::filter(genotype == "WT") %$% value,
  #        mito_tb_tidy %>% dplyr::filter(stage == "GV") %>% dplyr::filter(genotype == "Null") %$% value)
  # 
  # wilcox.test(mito_tb_tidy %>% dplyr::filter(stage == "MII") %>% dplyr::filter(genotype == "WT") %$% value,
  #        mito_tb_tidy %>% dplyr::filter(stage == "MII") %>% dplyr::filter(genotype == "Null") %$% value)
  
  # plot
  mito_plot <- 
    ggplot() +
    geom_jitter(data = mito_tb_tidy,
                mapping = aes(x = stage, y = value, color = genotype),
                size = 5, 
                position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4)) +
    stat_summary(data = mito_tb_tidy, 
                 mapping = aes(x = stage, y = value, color = genotype),
                 fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.4, 
                 position = position_dodge(width = 0.5)) +
    scale_fill_manual(values = c("#70ad47", "black")) +
    scale_color_manual(values = c("#70ad47", "black")) +
    scale_y_continuous(limits = c(0, max(mito_tb_tidy$value))) + 
    scale_x_discrete(drop = FALSE) +
    theme_bw() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.position = "none") +
    ggsave(filename = file.path(outpath, str_c("mitochondria_", mito_stat, ".png")), width = 7.5, height = 10)
  
  # return 
  return(mito_stat)
  
})


