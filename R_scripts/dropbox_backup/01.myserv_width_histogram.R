### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/MYSERV/FLI/golden_hamster.Siomi")

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

# genome path
genome_name <- basename(inpath)
genome_dir <- file.path(inpath, str_c("../../genomes/", genome_name))

# joined repeatMasker path
rmsk_path <- list.files(genome_dir, "rmsk\\..*\\.joined_rmsk_id\\.fa\\.out\\.gz", full.names = T)

# clean repeatMasker path rmsk.Siomi.20200701.clean.fa.out.gz
rmsk_clean_path <- list.files(genome_dir, "rmsk\\..*\\.clean\\.fa\\.out\\.gz", full.names = T)

######################################################## READ DATA
# read joined repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read clean repeatMasker
rmsk_clean <- readr::read_delim(rmsk_clean_path, delim = "\t")

######################################################## MAIN CODE
### get joined repeatMasker insertions
rmsk_filt <- 
  rmsk_tb %>% 
  dplyr::filter(str_detect(repName, "MYSERV"))

# get most common repName pairs
rmsk_repnames <- 
  rmsk_filt %$% 
  repName %>% 
  str_split(., "/") %>% 
  purrr::map(., function(repName) repName %>% unique(.) %>% sort(.) %>% str_c(., collapse = "/")) %>% 
  unlist(.) %>% 
  unique(.)

# get all MYSERVs
myserv_list <- 
  rmsk_clean$repName %>% 
  .[str_detect(., "MYSERV")] %>% 
  unique(.)

# get LTRs associated with MYSERVs
ltrs_list <- 
  rmsk_repnames[str_detect(rmsk_repnames, "/")] %>% 
  .[!str_detect(., "^IAP")] %>% 
  str_remove_all(., str_c(c(myserv_list, "/"), collapse = "|")) %>% 
  unique(.)
  
# get widths
rmsk_clean %<>%    
  dplyr::mutate(rmsk_id = as.character(rmsk_id), 
                width = end - start + 1)

### all IAPs
# filter 
rmsk_tb <- 
  rmsk_clean %>% 
  dplyr::filter(repName %in% c("MYSERV16_I", "MYSERV6-int",  "MYSERV-int", "RLTR31B2")) 

# plot histogram of widths
hist_width <-
  ggplot(data = rmsk_tb, aes(width, fill = repName)) +
  geom_histogram(binwidth = 10) +
  # scale_x_continuous(limits = c(100, 10000), breaks = seq(0, 10000, 100)) +
  facet_grid(rows = vars(repName), scales = "free_y") +
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
ggsave(filename = file.path(outpath, str_c("MYSERV", "RLTR31B2", "width_hist.facet.png", sep = ".")), plot = hist_width, width = 10, height = 5)
 