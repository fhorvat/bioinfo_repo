### INFO: 
### DATE: Sat Sep 28 15:06:52 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/datasets/2019_Sep/Analysis/piRNA_expression")

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

library(GenomicAlignments)
library(GenomicRanges)
library(GenomicFeatures)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# piRNA coordinates path
pirna_path <- file.path(genome_dir, "piRBase", "mmu.bed.gz")

# piRNA info path
pirna_info_path <- file.path(genome_dir, "piRBase", "piR_mmu.txt.gz")

# differentially expressed genes path
diffexp_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/datasets/2019_Sep/Analysis/expression/results/diffExp.DESeq2.genotype.all_biotype.significant_results.xlsx"

######################################################## READ DATA
# read piRNA coordinates
pirna_gr <- rtracklayer::import.bed(pirna_path)

# read piRNA info
pirna_info <- readr::read_delim(pirna_info_path, delim = "\t")
  
# read differentially expressed genes in DBL KO
diffexp_tb <- openxlsx::read.xlsx(diffexp_path, sheet = "DBL_vs_SOM")

######################################################## MAIN CODE
# filter piRNA info
pirna_info_tidy <- 
  pirna_info %>% 
  dplyr::filter(!is.na(aliases))

# filter piRNA clusters, reduce
pirna_gr_tidy <- 
  pirna_gr[mcols(pirna_gr)$name %in% pirna_info_tidy$name] %>% 
  reduce(., ignore.strand = F) 

# save as .bed for genomeBrowser
pirna_gr_tidy %>% 
  rtracklayer::export.bed(., file.path(pirna_path %>% str_replace("\\.bed\\.gz", ".annotated_piRNA.reduced.bed")))

# get coordinates of differentially expressed genes
diffexp_gr <- 
  diffexp_tb %>% 
  tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.)

# get overlaps between diff. expressed genes and piRNA clusters
diffexp_pirna_foverlaps <- findOverlaps(diffexp_gr, pirna_gr_tidy, ignore.strand = T)

# get piRNA clusters overlaping with differentially expressed genes
pirna_diffexp <- 
  pirna_gr_tidy[subjectHits(diffexp_pirna_foverlaps)] %>% 
  as_tibble(.) %>% 
  tidyr::unite(pirna_cluster, c("seqnames", "start", "end"), sep = " ")

# get differentially expressed genes which have annotated piRNA cluster inside
diffexp_pirna <- 
  diffexp_gr[queryHits(diffexp_pirna_foverlaps)] %>% 
  as_tibble(.) %>% 
  tidyr::unite(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::select(coordinates) %>% 
  dplyr::mutate(pirna_cluster = pirna_diffexp$pirna_cluster) %>% 
  dplyr::group_by(coordinates) %>% 
  dplyr::summarise(pirna_clusters = str_c(pirna_cluster, collapse = ", ")) %>% 
  dplyr::ungroup(.)

# add pirna clusters to original table
diffexp_tb_pirna_clusters <- 
  diffexp_tb %>% 
  as_tibble(.) %>% 
  dplyr::left_join(., diffexp_pirna, by = "coordinates") %>% 
  dplyr::filter(!is.na(pirna_clusters)) %T>%
  readr::write_csv(., path = file.path(outpath, "ensembl_genes_with_piRNA_clusters.DBL_vs_SOM.diffExp.csv"))
