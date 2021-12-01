### INFO: 
### DATE: Mon Aug 17 09:48:59 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.oocyte/map_clusters_to_Siomi")

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
library(BSgenome.Maur.UCSC.MesAur1)
library(GenomicAlignments)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# clusters table path
clusters_path <- file.path(inpath, "../expression/piRNA_clusters.oocyte_deduplexed.rpm_cutoff.10.20210518.csv")

# mapped clusters path
bam_path <- list.files(inpath, "*.bam$")

######################################################## READ DATA
# read cluster tables
clusters_tb <- readr::read_csv(clusters_path)

# get bam
bam <- GenomicAlignments::readGAlignmentsList(bam_path, use.names = T, param = ScanBamParam(what = "mapq"))
  
######################################################## MAIN CODE
# clean cluster table
clusters_tb_tidy <- 
  clusters_tb %>%
  dplyr::select(cluster_id, width_mesAur1 = width_full) %>% 
  dplyr::mutate(cluster_id = str_replace_all(cluster_id, " ", "_"))

# get table of alignments
# for each multimapping cluster get alignment with highest map quality
# calculate total length of split reads (splice = N + soft clip = S)
# filter by split length
bam_tb <- 
  bam %>% 
  unlist(.) %>% 
  as.data.frame(.) %>% 
  as_tibble(., rownames = "cluster_id") %>% 
  dplyr::mutate(cluster_id = str_remove(cluster_id, "\\.[0-9]+$")) %>% 
  dplyr::mutate(split_length = 
                  cigar %>% 
                  str_extract_all(., "[0-9]+N|[0-9]+S") %>% 
                  purrr::map(., function(x) str_remove_all(x, "N|S") %>% as.numeric(.) %>% sum) %>% 
                  unlist(.)) %>% 
  dplyr::group_by(cluster_id) %>% 
  dplyr::slice_max(order_by = mapq, n = 1, with_ties = TRUE) %>% 
  dplyr::slice_min(order_by = split_length, n = 1, with_ties = FALSE) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(split_length <= (qwidth / 2)) %>% 
  arrange(-qwidth)
  
# get GRanges for each cluster
clusters_gr <- GRanges(bam_tb)

# set names
names(clusters_gr) <- mcols(clusters_gr)$cluster_id

# save as table
readr::write_csv(as_tibble(clusters_gr), file.path(outpath, "piRNA_clusters.oocyte_deduplexed.rpm_cutoff.10.20210518.Siomi_coordinates.csv"))

# save as .bed
rtracklayer::export.bed(clusters_gr, file.path(outpath, "piRNA_clusters.oocyte_deduplexed.rpm_cutoff.10.20210518.Siomi_coordinates.bed"))



