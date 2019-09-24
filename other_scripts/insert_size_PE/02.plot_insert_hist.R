### INFO: 
### DATE: Tue Sep 24 13:35:46 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd(".")

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

# get paths of insert size histograms
hist_path <- list.files(inpath, pattern = "*\\.insertHis\\.txt", full.names = T)

######################################################## READ DATA
# read all histograms, join to one table
hist_tb <- purrr::map(hist_path, function(path){
  
  # read and clean
  readr::read_delim(path, delim = "\t", skip = 5) %>% 
    dplyr::mutate(sample_id = path %>% basename(.) %>% str_remove(., "\\.insertHis\\.txt"))
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  set_colnames(., c("insert_size", "count", "sample_id"))

# get mean values
mean_tb <- purrr::map(hist_path, function(path){
  
  # read and clean
  readr::read_lines(path, n_max = 1) %>% 
    str_extract(., "[0-9]+") %>% 
    tibble(mean = .) %>% 
    dplyr::mutate(mean = as.numeric(mean), 
                  sample_id = path %>% basename(.) %>% str_remove(., "\\.insertHis\\.txt"))
  
}) %>% 
  dplyr::bind_rows(.)

######################################################## MAIN CODE
# plot
plot_hist <- 
  ggplot() +
  geom_histogram(data = hist_tb, aes(x = insert_size, y = count, fill = sample_id, ), stat = "identity") +
  geom_vline(data = mean_tb, aes(xintercept = mean), linetype = 5) +
  geom_label(data = mean_tb, aes(x = mean, y = 0, label = mean)) +
  facet_grid(sample_id ~ .) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = - 0.2),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

# save plot
ggsave(filename = file.path(outpath, "insert_size.hist.pdf"), plot = plot_hist, 
       width = 20, height = 30)

