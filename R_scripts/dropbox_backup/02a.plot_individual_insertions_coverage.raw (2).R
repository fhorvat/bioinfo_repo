a### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/methylation/Bismark_nondirectional/individual_insertions")

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
library(rtracklayer)

library(ggrepel)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- file.path(getwd(), "..")

# set outpath
outpath <- getwd()

# mapped path
mapped_path <- inpath

# bigWigs path
bw_path <- list.files(inpath, ".*raw_coverage.bw$", full.names = T) 
bw_path <- bw_path[!str_detect(bw_path, "PE_s")]

# methylation data path
meth_path <- list.files(inpath, ".*\\.bismark\\.cov\\.gz", full.names = T)

######################################################## READ DATA
# read coverage data
coverage_tb_full <- purrr::map(bw_path, function(path){
  
  # read bw
  bw_gr <- rtracklayer::import.bw(path, as = "RleList")
  
  # get coverage
  coverage_list <- 
    bw_gr %>% 
    as(., "IntegerList") %>% 
    unlist(.) %>% 
    split(., names(.)) 
  
  # set to table
  coverage_tb <- purrr::map(names(coverage_list), function(seqname){
    
    # create table
    cov_tb <- 
      coverage_list[[seqname]] %>% 
      tibble(coverage = .) %>% 
      dplyr::mutate(repName = seqname,
                    sample_id = basename(path), 
                    pos = 1:nrow(.)) %>% 
      dplyr::select(repName, sample_id, pos, coverage)
    
  }) %>% 
    dplyr::bind_rows(.)
  
  # return
  return(coverage_tb)
  
}) %>% 
  dplyr::bind_rows(.)


# read methylation data
meth_tb_full <- purrr::map(meth_path, function(path){
  
  # read table
  meth_tb <- 
    path %>% 
    readr::read_delim(., delim = "\t", col_names = c("repName", "start", "end", "percentage", "methylated", "unmethylated")) %>% 
    dplyr::mutate(sample_id = path %>% basename(.) %>% stringr::str_remove(., "\\.converted_bismark_bt2\\.bismark\\.cov\\.gz"))
  
}) %>% 
  dplyr::bind_rows(.)

######################################################## MAIN CODE
### get data for plot 
# filter table
coverage_tb <- 
  coverage_tb_full %>% 
  dplyr::mutate(sample_id = str_remove_all(sample_id, "\\.converted_bismark_bt2.raw_coverage.bw"))

# filter methylation table - don't take in account CpG sites with less than 30 reads
meth_tb <- 
  meth_tb_full %>% 
  dplyr::filter((methylated + unmethylated) >= 30)
  
# polygon total coverage
pg <- tibble(px = coverage_tb$pos, 
             py = coverage_tb$coverage, 
             sample_id = coverage_tb$sample_id, 
             repName = coverage_tb$repName) 

# polygon methylation
pg_meth <- tibble(px = meth_tb$start, 
                  percentage = meth_tb$percentage, 
                  sample_id = meth_tb$sample_id, 
                  repName = meth_tb$repName)

# loop through individual insertions
for(insertion in unique(pg$repName)){
  
  # filter coverage table
  pg_filt <- 
    pg %>% 
    dplyr::filter(repName == insertion)
  
  # filter methylation table
  pg_meth_tb <- 
    pg_meth %>% 
    dplyr::filter(repName == insertion) %>% 
    dplyr::mutate(py = (percentage * 0.01 * max(pg_filt$py)), 
                  percentage = round(percentage, 2))
  
  # create plot
  baseplot <- 
    ggplot() +
    geom_rect(data = pg_filt, aes(xmin = px, xmax = px + 1, ymin = 0, ymax = py), fill = "black") +
    geom_rect(data = pg_meth_tb, aes(xmin = px, xmax = px + 1, ymin = 0, ymax = py), color = "red", alpha = 1) +
    geom_text_repel(data = pg_meth_tb, aes(x = px, y = py, label = percentage), force = 2, angle = 90, direction = "x", segment.color = NA) +
    facet_grid(cols = vars(sample_id), 
               scales = "free_y") +
    ggtitle(str_replace_all(insertion, ":|-", " ")) + 
    ylab("") +
    xlab("") +
    # scale_x_continuous(limits = c(0, 6500)) +
    # scale_y_continuous(limits = c(0, 90)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          # plot.title = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save plot
  ggsave(plot = baseplot,
         filename = file.path(outpath, "test.png"),
         width = 15, height = 7.5)
  
  # # save plot
  # ggsave(plot = baseplot, 
  #        filename = file.path(outpath, str_c("raw_coverage", insertion, "bw", "jpg", sep = ".")), 
  #        width = 8, height = 2)
  
}
