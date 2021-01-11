### INFO: 
### DATE: Mon Aug 17 09:48:59 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/map_clusters_to_Siomi")

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
clusters_path <- file.path(inpath, "MesAur1.1k_pachytene_clusters.200730.recalculated.fused.PS.xlsx")

# mapped clusters path
bam_path <- file.path(inpath, "MesAur1.1k_pachytene_clusters.200730.recalculated.fused.PS.bam")

######################################################## READ DATA
# read cluster tables
clusters_tb <- 
  openxlsx::read.xlsx(clusters_path) %>% 
  as_tibble(., .name_repair = "unique")

# get bam
bam <- GenomicAlignments::readGAlignmentsList(bam_path, use.names = T)
  
######################################################## MAIN CODE
# get alignment for each cluster
clusters_gr <- 
  bam %>% 
  grglist(.) %>% 
  range(.) %>% 
  unlist(.)

# get names as coordinates
mcols(clusters_gr)$coordinates <- names(clusters_gr)

# get in table
cluster_tb_Siomi <- 
  clusters_gr %>% 
  as_tibble(.) %>% 
  dplyr::group_by(coordinates) %>% 
  dplyr::filter(width == max(width)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::unite(coordinates.Siomi, seqnames, start, end, sep = " ") %>% 
  dplyr::select(coordinates, coordinates.Siomi) %>% 
  dplyr::mutate(coordinates = str_replace_all(coordinates, "_", " "))

# join with original table
clusters_tb_joined <- 
  clusters_tb %>% 
  dplyr::left_join(., cluster_tb_Siomi, by = "coordinates") %>% 
  dplyr::select(coordinates, coordinates.Siomi, everything()) 

# save
# openxlsx::write.xlsx(clusters_tb_joined, "MesAur1.1k_pachytene_clusters.200730.recalculated.fused.Siomi_coordinates.FH.xlsx")

### get clusters on chromosome X
clusters_x <- 
  clusters_tb_joined %>% 
  dplyr::filter(str_detect(coordinates.Siomi, "HiC_scaffold_22"))