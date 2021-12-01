### INFO: 
### DATE: Mon Aug 17 09:48:59 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.oocyte/map_clusters_to_Siomi/LINE1_FLI_expression")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# clusters table path
clusters_path <- file.path(inpath, "../piRNA_clusters.oocyte_deduplexed.rpm_cutoff.10.20210518.Siomi_coordinates.csv")

# library size path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.deduplicated.smallRNAseq/Data/Mapped/STAR_Siomi.multimappers"
library_size_path <- file.path(mapped_path, "4_library_size/library_sizes.txt")

# mapped clusters path
bam_path <- file.path(inpath, "03_bam_clusters_subset")
bam_path <- list.files(bam_path, "*.bam$", full.names = T)

######################################################## READ DATA
# read cluster tables
clusters_tb <- readr::read_csv(clusters_path)

# read library size table
library_size_tb <- readr::read_delim(library_size_path, delim = "\t", col_names = c("sample_id", "library_size"))

# get bam
bam_list <- purrr::map(bam_path, function(path){
  
  # reads with NH tag
  GenomicAlignments::readGAlignmentsList(path, use.names = T, param = ScanBamParam(tag = "NH")) %>% 
    unlist(.)
  
}) %>% 
  set_names(basename(bam_path) %>% str_remove(., "\\.bam"))

######################################################## MAIN CODE
# get clusters GRanges
cluster_gr <- GRanges(clusters_tb)

# summarize overlaps between reads and clusters
sum_overlaps <- purrr::map(names(bam_list), function(sample_id){
  
  # get one bam
  bam <- bam_list[[sample_id]]
  
  # find overlaps
  overlaps <- findOverlaps(cluster_gr, bam, ignore.strand = T)
  
  # for each cluster get count of reads, both full reads and fractions for multimappers
  count_tb <- 
    tibble(cluster_id = cluster_gr[queryHits(overlaps)]$cluster_id, 
           count_full = rep(1, length(subjectHits(overlaps))), 
           count_fraction = (1 / mcols((bam[subjectHits(overlaps)]))$NH)) %>% 
    dplyr::group_by(cluster_id) %>% 
    dplyr::summarise(count_full = sum(count_full), 
                     count_fraction = sum(count_fraction)) %>% 
    dplyr::mutate(sample_id = sample_id)
  
  # return 
  return(count_tb)
  
}) %>% 
  dplyr::bind_rows(.)

# get RPM
rpm_tb <-
  sum_overlaps %>% 
  dplyr::left_join(., library_size_tb, by = "sample_id") %>% 
  dplyr::mutate(library_size = (library_size / 1e6), 
                rpm_full = count_full / library_size, 
                rpm_fraction = count_fraction / library_size) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "WT|KO")) %>% 
  dplyr::group_by(cluster_id, genotype) %>% 
  dplyr::summarise(rpm_full = mean(rpm_full), 
                   rpm_fraction = mean(rpm_fraction)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::pivot_wider(id_cols = cluster_id, names_from = genotype, values_from = c(rpm_full, rpm_fraction), values_fill = 0, names_sep = ".") %>% 
  dplyr::arrange(-rpm_full.KO) %>% 
  dplyr::left_join(., clusters_tb %>% 
                     dplyr::select(cluster_id, seqnames:end) %>% 
                     tidyr::unite(col = coordinates.Siomi, seqnames, start, end, sep = " ", remove = T), 
                   by = "cluster_id") %>% 
  dplyr::mutate(coordinates = str_remove(cluster_id, "piwil[1,3]_") %>% str_replace_all(., "_", " ")) %>% 
  dplyr::select(cluster_id, coordinates.mesAur1 = coordinates, coordinates.Siomi, everything())

# write
readr::write_csv(rpm_tb, file.path(outpath, "piRNA_clusters.oocyte_deduplexed.Siomi_LINE1_FLI_reads.RPM.csv"))



