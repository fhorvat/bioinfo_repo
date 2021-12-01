### INFO: 
### DATE: Sun Sep 08 20:48:52 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/ATP_concentration")

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

# ATP table path
atp_path <- file.path(inpath, "final_atp_assay_results_SG.190908PS.csv")

######################################################## READ DATA
# read ATP table
atp_tb <- readr::read_csv(atp_path)

######################################################## MAIN CODE
# gather + separate table
atp_tb_tidy <- 
  atp_tb %>% 
  tidyr::gather(key = "timepoint", value = "value") %>% 
  tidyr::separate(timepoint, into = c("genotype", "timepoint"), sep = " ") %>% 
  dplyr::filter(!is.na(value)) %>% 
  dplyr::mutate(genotype = factor(genotype, levels = c("WT", "KO")))

wilcox.test(atp_tb_tidy %>% dplyr::filter(timepoint == "0h") %>% dplyr::filter(genotype == "WT") %$% value,
            atp_tb_tidy %>% dplyr::filter(timepoint == "0h") %>% dplyr::filter(genotype == "KO") %$% value)

wilcox.test(atp_tb_tidy %>% dplyr::filter(timepoint == "20h") %>% dplyr::filter(genotype == "WT") %$% value,
            atp_tb_tidy %>% dplyr::filter(timepoint == "20h") %>% dplyr::filter(genotype == "KO") %$% value)

# plot
atp_plot <- 
  ggplot() +
  geom_jitter(data = atp_tb_tidy,
              mapping = aes(x = timepoint, y = value, color = genotype),
              size = 5, 
              position = position_jitterdodge(dodge.width = 1, jitter.width = 0.6)) +
  stat_summary(data = atp_tb_tidy, 
               mapping = aes(x = timepoint, y = value, color = genotype),
               fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.8, 
               position = position_dodge(width = 1)) +
  scale_fill_manual(values = c("#70ad47", "black")) +
  scale_color_manual(values = c("#70ad47", "black")) +
  scale_y_continuous(limits = c(0, max(atp_tb_tidy$value))) + 
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none") +
  ggsave(filename = file.path(outpath, "atp_concentration.jitter.2.png"), width = 7.5, height = 10)



