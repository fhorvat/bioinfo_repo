### INFO: 
### DATE: Mon Aug 31 14:57:45 2020
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

library(openxlsx)
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

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# arguments from command line
clusters_path <- args$clusters_path
bam_path <- args$bam_path

# bam path
bam_name <- basename(bam_path) %>% str_remove(., "\\.bam")

######################################################## READ DATA
# read clusters tb
clusters_tb <- 
  openxlsx::read.xlsx(clusters_path) %>% 
  as_tibble(., .name_repair = "unique")

######################################################## MAIN CODE
### tidy and prepare data
# clean clusters table, get GRanges
clusters_gr <- 
  clusters_tb %>% 
  dplyr::select(coordinates) %>% 
  dplyr::filter(!is.na(coordinates)) %>% 
  tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ", remove = F) %>% 
  GenomicRanges::GRanges(.)

# read bam
bam_gr <- 
  GenomicAlignments::readGAlignmentsList(file = bam_path, 
                                         use.names = TRUE, 
                                         param = ScanBamParam(which = clusters_gr, 
                                                              flag = scanBamFlag(isMinusStrand = NA))) %>% 
  unlist(.)

### get coverage of all clusters
# create GRangesList
clusters_grlist <- split(clusters_gr, clusters_gr$coordinates)

# subset coverage of both strands for each cluster - separate for + and - strands
coverage_tb <- lapply(clusters_grlist, function(cluster_gr){
  
  # get length of chromosome on which is feature located
  seq_length <- seqlengths(bam_gr)[as.character(seqnames(cluster_gr))]
  
  # subset bam for reads overlapping cluster
  cluster_gr_subset <- subsetByOverlaps(bam_gr, cluster_gr, ignore.strand = T)
  
  # set seqnames to only cluster chromosome/scaffold
  seqlevels(cluster_gr_subset) <- as.character(unique(seqnames(cluster_gr_subset)))

  #  get coverage on separate strands
  coverage_stranded <- purrr::map(c("+", "-"), function(strand_sub){
    
    # get coverage
    cluster_coverage <- 
      cluster_gr_subset[strand(cluster_gr_subset) == strand_sub] %>% 
      coverage(.) %>% 
      as(., "IntegerList") %>% 
      unlist(.) %>% 
      unname(.)
    
    if(length(cluster_coverage) == 0){
      coverage_df <- tibble(pos = 1:seq_length, 
                            coverage = 0)
    }else{
      coverage_df <- tibble(pos = 1:seq_length, 
                            coverage = cluster_coverage)
    }
    
    # set position to 0, add strand info
    coverage_df %<>% 
      dplyr::filter((pos >= start(cluster_gr)) & (pos <= end(cluster_gr))) %>% 
      dplyr::mutate(pos = 1:nrow(.), 
                    strand = strand_sub)
    
    # return
    return(coverage_df)
    
  }) %>% 
    dplyr::bind_rows(.) %>% 
    dplyr::mutate(cluster = mcols(cluster_gr)$coordinates)
  
  # return 
  return(coverage_stranded)
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(sample_id = bam_name) %>%
  dplyr::select(sample_id, cluster, pos, coverage, strand)

# save 
readr::write_csv(coverage_tb, file.path(outpath, str_c(bam_name, "MesAur1.1k_pachytene_clusters.200730.coverage.csv", sep = ".")))
