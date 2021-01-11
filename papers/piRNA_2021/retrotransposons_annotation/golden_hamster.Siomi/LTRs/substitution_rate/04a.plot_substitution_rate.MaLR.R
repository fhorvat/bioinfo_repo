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
# set classes
comboClasses <- list("MLT1O" = c("MLT1O"), 
                     "MLT1N" = c("MLT1N2"), 
                     "MLT1M" = c("MLT1M"), 
                     "MLT1L" = c("MLT1L"),
                     "MLT1K" = c("MLT1K"), 
                     "MLT1J" = c("MLT1J", "MLT1J1", "MLT1J2"), 
                     "MLT1I" = c("MLT1I"), 
                     "MLT1H" = c("MLT1H", "MLT1H1", "MLT1H2"), 
                     "MLT1G" = c("MLT1G", "MLT1G1", "MLT1G3"),
                     "MLT1F" = c("MLT1F", "MLT1F1", "MLT1F2"), 
                     "MLT1E" = c("MLT1E", "MLT1E1", "MLT1E1A", "MLT1E2", "MLT1E3"), 
                     "MLT1D" = c("MLT1D"), 
                     "MLT1C" = c("MLT1C"), 
                     "MLT1B" = c("MLT1B"),
                     "MLT1A" = c("MLT1A", "MLT1A0", "MLT1A1"), 
                     "MLT2F" = c("MLT2F"),
                     "MLT2E" = c("MLT2E"),
                     "MLT2D" = c("MLT2D"),
                     "MLT2C" = c("MLT2C1", "MLT2C2"), 
                     "MLT2B" = c("MLT2B1", "MLT2B2", "MLT2B3", "MLT2B4", "MLT2B5"), 
                     "ORR1G" = c("ORR1G"),
                     "ORR1F" = c("ORR1F"),
                     "ORR1E" = c("ORR1E"),
                     "ORR1D" = c("ORR1D1", "ORR1D2"),
                     "ORR1C" = c("ORR1C1", "ORR1C2"),
                     "ORR1B" = c("ORR1B1", "ORR1B2"),
                     "ORR1A" = c("ORR1A0", "ORR1A1", "ORR1A2", "ORR1A3", "ORR1A4"),
                     "MTE"   = c("MTEa", "MTEb", "MTE2a", "MTE2b"), 
                     "MTD"   = c("MTD"),
                     "MTC"   = c("MTC"),
                     "MTB"   = c("MTB", "MTB_Mm"),
                     "MTA"   = c("MTA_Mm"),
                     "MT2C"  = c("MT2C_Mm"),
                     "MT2B"  = c("MT2B", "MT2B1", "MT2B2"),
                     "MT2A"  = c("MT2A"),
                     "MT2"   = c("MT2_Mm"))



# list to tibble
classes_tb <- purrr::map(names(comboClasses), function(ltr_name){
  
  # get one LTR type in tibble
  comboClasses[[ltr_name]] %>% 
    tibble(repName = ., repSubfamily = ltr_name)
  
}) %>%
  dplyr::bind_rows(.)

# remove ORR1A and MTA from table
classes_tb %<>% 
  dplyr::filter(!(repSubfamily %in% c("ORR1A", "MTA")))


### plot
# prepare table for plot - add category and repFamily, resample to include only 200 per each subFamily
sub_rate_tb <- 
  sub_rate_tb_all %<>% 
  dplyr::right_join(., classes_tb, by = "repName") %>%
  dplyr::mutate(repSubfamily = factor(repSubfamily, levels = unique(classes_tb$repSubfamily))) %>% 
  dplyr::group_by(repSubfamily) %>% 
  sample_n(ifelse(n() >= 200, 200, n())) %>% 
  dplyr::ungroup(.)

# data statistics
sub_rate_stats <- 
  sub_rate_tb %>% 
  group_by(repSubfamily) %>% 
  summarize(median = median(sub_rate), 
            sd = sd(sub_rate), 
            length = length(sub_rate)) %>% 
  dplyr::mutate(SEM = sd / sqrt(length))

# labels for plot
labels_tb <- 
  sub_rate_tb %>% 
  dplyr::group_by(repSubfamily) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::filter(count > 0)

# plot
sub_rate_boxplot <- 
  ggplot() +
  stat_boxplot(data = sub_rate_tb, aes(x = repSubfamily, y = sub_rate), geom = "errorbar", size = 1.5) +
  geom_boxplot(data = sub_rate_tb, aes(x = repSubfamily, y = sub_rate), outlier.colour = NULL, outlier.shape = NA, size = 1.5) +
  geom_jitter(data = sub_rate_tb, aes(x = repSubfamily, y = sub_rate), alpha = 0.4, 
              colour = "black", size = 1.5, width = 0.1, height = 0, show.legend = F) +
  scale_x_discrete(labels = str_c(labels_tb$repSubfamily, 
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
       filename = file.path(outpath, str_c("LTRs.random_200_per_subFamily.200730", "MaLR", "boxplot.png", sep = ".")), 
       width = 20, 
       height = 10)

