### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/IAP")

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

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# # joined repeatMasker path
# rmsk_path <- file.path(genome_dir, "rmsk.Siomi.20200701.joined_rmsk_id.fa.out.gz")

# clean repeatMasker path
rmsk_clean_path <- file.path(genome_dir, "rmsk.Siomi.20200701.clean.fa.out.gz")

######################################################## READ DATA
# # read joined repeatMasker
# rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read clean repeatMasker
rmsk_clean <- readr::read_delim(rmsk_clean_path, delim = "\t")

######################################################## MAIN CODE
# get widths
rmsk_clean %<>%    
  dplyr::mutate(rmsk_id = as.character(rmsk_id), 
                width = end - start + 1)

### all IAPs
# filter 
rmsk_tb <- 
  rmsk_clean %>% 
  dplyr::filter(str_detect(repName, "^IAP")) 

# plot histogram of widths
hist_width <-
  ggplot(data = rmsk_tb, aes(width), color = "black") +
  geom_histogram(binwidth = 10) +
  # scale_x_continuous(limits = c(100, 10000), breaks = seq(0, 10000, 100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

# save
ggsave(filename = file.path(outpath, str_c("01", "IAP", "all", "width_hist.facet.png", sep = ".")), plot = hist_width, width = 10, height = 5)


### LTR IAPs
# filter 
rmsk_tb <- 
  rmsk_clean %>% 
  dplyr::filter(str_detect(repName, "^IAP")) %>% 
  dplyr::filter(!(str_detect(repName, "-int|_I$")))

# plot histogram of widths
hist_width <-
  ggplot(data = rmsk_tb, aes(width), color = "black") +
  geom_histogram(binwidth = 10) +
  # scale_x_continuous(limits = c(100, 10000), breaks = seq(0, 10000, 100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

# save
ggsave(filename = file.path(outpath, str_c("02", "IAP", "LTRs", "width_hist.facet.png", sep = ".")), plot = hist_width, width = 10, height = 5)


### Int IAPs
# filter 
rmsk_tb <- 
  rmsk_clean %>% 
  dplyr::filter(str_detect(repName, "^IAP")) %>% 
  dplyr::filter((str_detect(repName, "-int|_I$")))

# plot histogram of widths
hist_width <-
  ggplot(data = rmsk_tb, aes(width), color = "black") +
  geom_histogram(binwidth = 10) +
  # scale_x_continuous(limits = c(100, 10000), breaks = seq(0, 10000, 100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

# save
ggsave(filename = file.path(outpath, str_c("03", "IAP", "internal", "width_hist.facet.png", sep = ".")), plot = hist_width, width = 10, height = 5)


### selection IAPLTR3+IAPLTR4 LTR+int
# filter 
rmsk_tb <- 
  rmsk_clean %>% 
  dplyr::filter(str_detect(repName, "^IAP")) %>% 
  dplyr::filter(repName %in% c("IAPLTR3", "IAPLTR3-int", "IAPLTR4", "IAPLTR4_I"))

# plot histogram of widths
hist_width <-
  ggplot(data = rmsk_tb, aes(width), color = "black") +
  facet_wrap(vars(repName), nrow = 4, ncol = 1) +
  geom_histogram(binwidth = 10) +
  # scale_x_continuous(limits = c(100, 10000), breaks = seq(0, 10000, 100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

# save
ggsave(filename = file.path(outpath, str_c("04", "IAP", "IAPLTR3_4.LTRs_ints", "width_hist.facet.png", sep = ".")), plot = hist_width, width = 10, height = 10)


### intact guys
# read intact table
fli_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/IAP"
fli_path <- file.path(fli_path, "IAP.FLI_elements.csv")
rmsk_tb <- readr::read_csv(fli_path)

# plot histogram of widths
hist_width <-
  ggplot(data = rmsk_tb, aes(width), color = "black") +
  geom_histogram(binwidth = 5) +
  # scale_x_continuous(limits = c(0, max(rmsk_tb$width)) + 100) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

# save
ggsave(filename = file.path(outpath, str_c("05", "IAP", "intact_guys", "width_hist.facet.png", sep = ".")), plot = hist_width, width = 10, height = 10)
