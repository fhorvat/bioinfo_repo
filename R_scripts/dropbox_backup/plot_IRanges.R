### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/PhD/algorithms_and_programming/2020_10_26/GRanges/lectures")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS
IRangesMethodToTable <- function(IRangesMethod, ...){
  
  # apply method to IRanges
  irange_result <- IRangesMethod(...)
  
  # create list of IRanges
  irange_list <- list(..., irange_result)
  
  # loop through the list: transform to table, add unique id, bind to one table
  irange_tb <- lapply(1:length(irange_list), function(index){
    
    if(class(irange_list[[index]]) %in% c("IRanges", "GRanges")){
      
      irange_list[[index]] %>% 
        as.data.frame(.) %>% 
        as_tibble(.) %>% 
        dplyr::mutate(range_id = index)
      
    }
    
  }) %>% 
    dplyr::bind_rows(.) %>% 
    
    # replace last added range with result
    dplyr::mutate(range_id = ifelse(test = (range_id == max(range_id)), 
                                    yes = "result", 
                                    no = paste("range", range_id, sep = ".")))
  
  # return table
  return(irange_tb)
  
}

plotIRanges <- function(irange_tb){
  
  # calculate max. x-limit
  xlim_max <- max(irange_tb$end) + 1
  
  # add range ID if it doesn't exist
  if(!("range_id" %in% names(irange_tb))){
    
    irange_tb %<>% 
      dplyr::mutate(range_id = "range.1")
    
  }
  
  # add bins to table
  irange_tb_bins <- purrr::map(unique(irange_tb$range_id), function(id){
    
    # get table
    bins_tb <- 
      irange_tb %>% 
      dplyr::filter(range_id == id)
    
    # get IRanges object, disjoint bins
    bins <- 
      IRanges(start = bins_tb$start, end = bins_tb$end + 1) %>% 
      disjointBins(.)
    
    # add bins to table
    bins_tb %<>% 
      dplyr::mutate(range_bin = bins)
    
    # return
    return(bins_tb)
    
  }) %>% 
    dplyr::bind_rows(.)
  
  # add range ID group, modify bin value by group
  irange_tb_bins_plot <- 
    irange_tb_bins %>% 
    dplyr::mutate(range_id = factor(range_id, levels = rev(unique(range_id)))) %>% 
    dplyr::group_by(range_id) %>% 
    dplyr::mutate(range_group = dplyr::cur_group_id()) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::mutate(range_bin = ((range_group * length(unique(.$range_bin)) + 1) + range_bin)) %>% 
    dplyr::mutate(range_id = as.character(range_id)) %>% 
    dplyr::mutate(range_id = replace(range_id, range_id == "result", NA))
  
  # plot
  ggplot() +
    geom_rect(data = irange_tb_bins_plot, aes(xmin = start, xmax = end,
                                              ymin = range_bin, ymax = range_bin + 0.9,
                                              fill = range_id)) +
    scale_fill_grey(start = 0.8, end = 0.3, na.value = "red3") +
    scale_x_continuous(breaks = 0:xlim_max) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = "none")
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

######################################################## READ DATA

######################################################## MAIN CODE
# create ranges
range_1 <- IRanges(start = c(4, 18, 1, 2), 
                   width = c(3, 4, 3, 5))

range_2 <- IRanges(start = c(8, 3, 16, 3), 
                   width = c(5, 10, 4, 7))

IRangesMethodToTable(IRanges::punion, range_1, range_2, fill.gap = T) %>% plotIRanges(.)
