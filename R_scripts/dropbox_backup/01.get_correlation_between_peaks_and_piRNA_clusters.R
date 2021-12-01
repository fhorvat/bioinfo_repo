### INFO: 
### DATE: Thu Jul 23 13:37:16 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/compare_with_chipseq/Walker_2015_EpigenetChromatin_GSE61613")

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

library(GenomicDistributions)
library(GenomicRanges)
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# peak caller path
peaks_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/ChIPseq/Walker_2015_EpigenetChromatin_GSE61613/Analysis/peak_call"

# peak caller narrowPeak files
peaks_narrow_path <- 
  list.files(peaks_path, pattern = "\\.narrowPeak", full.names = T, recursive = T) %>% 
  .[!str_detect(., "UCSC")]

# cluster coordinates path
cluster_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis"
cluster_tb_path <- list.files(cluster_path, pattern = "\\.xlsx", full.names = T)

######################################################## READ DATA
# read narrow peak files
peaks_narrow_grlist <- purrr::map(peaks_narrow_path, function(path){
  
  # read table,
  peaks_gr <- 
    path %>% 
    readr::read_delim(., delim = "\t", 
                      col_names = c("seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")) %>% 
    GRanges(.) %>% 
    sort(.)
  
}) %>%
  set_names(., str_extract(peaks_narrow_path, "ChIP_peaks|affinity_peaks")) %>% 
  GRangesList(.)

# read clusters 
clusters_tb <- openxlsx::read.xlsx(cluster_tb_path, startRow = 2)

######################################################## MAIN CODE
# get clusters as GRanges
clusters_gr <- 
  clusters_tb %>% 
  as_tibble(.) %>% 
  dplyr::select(coordinates = refined.mouse.piRNA.cluster) %>% 
  tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = ":|-") %>% 
  dplyr::mutate_at(vars("start", "end"), ~(str_replace_all(., ",", "") %>% as.numeric(.))) %>% 
  dplyr::filter_all(all_vars(!is.na(.))) %>% 
  GRanges(.)

# calculate distance from peaks to clusters
peaks_clusters_dist <- calcFeatureDist(peaks_narrow_grlist, clusters_gr)

# print heatmap
png(filename = file.path(outpath, "peaks_to_piRNA_cluster.distance.hist.png"), width = 10, height = 10, units = "in", res = 300)
plotFeatureDist(peaks_clusters_dist, featureName = "piRNA clusters", tile = FALSE, nbin = 20, numbers = T, size = 200000) 
dev.off()






