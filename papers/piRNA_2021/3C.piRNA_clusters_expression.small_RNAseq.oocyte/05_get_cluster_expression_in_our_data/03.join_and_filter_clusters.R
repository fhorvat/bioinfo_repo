### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.oocyte/oocyte_expression")

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
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# clusters path
clusters_list <- list.files(inpath, ".*RPKM.deduplexed.our_data.csv$", full.names = T, recursive = T)

######################################################## READ DATA
# 10 RPM minimum in one or the other IP and 
# 10 RPM minimum for 18-20 & 24-32 reads combined in WT average.

# read and filter clusters
clusters_filtered_list <- purrr::map(clusters_list, function(path){
  
  # read cluster table
  cluster_tb <- 
    readr::read_csv(path) %>% 
    dplyr::filter_at(vars(starts_with("rpm_PIWIL")), dplyr::any_vars(. >= 10)) %>% 
    dplyr::filter((rpm_Mov10l1_WT.18to20nt + rpm_Mov10l1_WT.24to27nt + rpm_Mov10l1_WT.28to31nt) >= 10) %>% 
    dplyr::mutate(cluster_defined = str_extract(path, "piwil1|piwil3"), 
                  cluster_id = str_c(cluster_defined, coordinates, sep = "_")) %>% 
    dplyr::select(cluster_id, everything())
  
  # return
  return(cluster_tb)
  
}) %>% 
  set_names(str_extract(clusters_list, "piwil1|piwil3"))

### get overlap between PIWIL1 and PIWIL3 defined clusters
# get GRanges
clusters_gr_list <- purrr::map(clusters_filtered_list, function(cluster_tb){
  
  # get table, transform to GRanges
  clusters_gr <- 
    cluster_tb %>% 
    dplyr::select(cluster_id, coordinates, width_full) %>% 
    tidyr::separate(col = coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
    GRanges(.)
  
  # return 
  return(clusters_gr)
  
})

# overlap
piwil1_gr <- clusters_gr_list[["piwil1"]]
piwil3_gr <- clusters_gr_list[["piwil3"]]
overlaps <- findOverlaps(piwil1_gr, piwil3_gr, ignore.strand = T)

# get clusters which overlap and merge them
clusters_piwil_both <- 
  c(piwil1_gr[queryHits(overlaps)], piwil3_gr[subjectHits(overlaps)]) %>% 
  GenomicRanges::reduce(., ignore.strand = T)
mcols(clusters_piwil_both)$cluster_in <- c("PIWIL1/PIWIL3")

# get clusters which don't overlap
clusters_piwil1 <- piwil1_gr[setdiff(1:length(piwil1_gr), queryHits(overlaps))]
mcols(clusters_piwil1) <- NULL
mcols(clusters_piwil1)$cluster_in  <- "PIWIL1"

clusters_piwil3 <- piwil3_gr[setdiff(1:length(piwil3_gr), subjectHits(overlaps))]
mcols(clusters_piwil3) <- NULL
mcols(clusters_piwil3)$cluster_in  <- "PIWIL3"

# join GRanges
clusters_piwil_all <- c(clusters_piwil_both, clusters_piwil1, clusters_piwil3)

# set names
names(clusters_piwil_all) <- str_c(seqnames(clusters_piwil_all), ":", start(clusters_piwil_all), "-", end(clusters_piwil_all))

# export as .gtf
export.gff3(object = clusters_piwil_all, 
            con = file.path(outpath, 
                            "piRNA_clusters.oocyte_deduplexed.rpm_cutoff.10.20210518.gff"))

# export as .bed
export.bed(object = clusters_piwil_all, 
           con = file.path(outpath, 
                           "piRNA_clusters.oocyte_deduplexed.rpm_cutoff.10.20210518.bed"))

# create table
clusters_piwil_all_tb <- 
  clusters_piwil_all %>% 
  as_tibble(.) %>% 
  dplyr::arrange(-width)

# save
readr::write_csv(clusters_piwil_all_tb, file.path(outpath, 
                                                  "piRNA_clusters.oocyte_deduplexed.rpm_cutoff.10.20210518.csv"))



