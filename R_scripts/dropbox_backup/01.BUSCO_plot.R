### INFO: 
### DATE: Mon Sep 02 09:57:40 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/Analysis")

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
suppressAll <- function(x) suppressMessages(suppressWarnings(x))

######################################################## PATH VARIABLES
# set inpath 
inpath <- file.path(getwd(), "../transcriptome_assemblies")

# set outpath
outpath <- getwd()

# find Busco reports
busco_path <- list.files(inpath, pattern = "short_summary_s_.*\\.GG\\.Trinity\\.busco\\.mammalia\\.txt", recursive = T, full.names = T)

######################################################## READ DATA
# read and clean BUSCO report
busco_tb <- purrr::map(busco_path, function(path){
  
  tb <- 
    suppressAll(readr::read_delim(path, delim = "\t", col_names = c("tmp", "count", "category"))) %>% 
    dplyr::select(category, count) %>% 
    dplyr::mutate(category = ifelse(str_detect(count, "C:[0-9]+\\.[0-9]+\\%\\[.*"), count, category)) %>% 
    dplyr::filter_all(all_vars(!is.na(.))) %>% 
    dplyr::mutate(category = ifelse(str_detect(category, "C:[0-9]+\\.[0-9]+\\%\\[.*"), "Percentage of BUSCOs found", category), 
                  count = ifelse(str_detect(count, "C:[0-9]+\\.[0-9]+\\%\\[.*"), count %>% str_remove_all(., "C:|\\%\\[.*"), count), 
                  count = as.numeric(count), 
                  animal = path %>% str_remove_all(., "/common/.*/transcriptome_assemblies/|trinity_|\\.GG") %>% str_remove(., "\\..*"))
  
}) %>% 
  dplyr::bind_rows(.)

######################################################## MAIN CODE
# tidy table for plot
plot_tb <- 
  busco_tb %>% 
  dplyr::filter(category == "Percentage of BUSCOs found") %>% 
  dplyr::select(-category) %>% 
  dplyr::mutate(animal = factor(animal, 
                                levels = c("mouse_B6", "mouse_PWD", "mouse_cast", "rat", 
                                           "golden_hamster", "chinese_hamster", "rabbit", 
                                           "pig", "cow"),
                                labels = c("Mouse CL57/BL6", "Mouse PWD", "Mouse CAST/EiJ", "Rat", 
                                           "Golden hamster", "Chinese hamster", "Rabbit", 
                                           "Pig", "Cow"))) %>% 
  dplyr::arrange(animal)

# plot barplot
barplot <- 
  ggplot(data = plot_tb, aes(x = animal, y = count, fill = animal)) +
  geom_bar(stat = "identity") +
  # geom_hline(yintercept = 100) +
  scale_y_continuous(limits = c(0, 100)) +
  ylab("Percentage of BUSCOs found") +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 15)) +
  # theme(panel.border = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15, vjust = 0.5), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none")

# save plot
ggsave(plot = barplot, filename = file.path(outpath, "busco_percentage.mammalian.png"), width = 12, height = 10)

