### INFO: 
### DATE: Thu May 02 16:22:37 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/other/lab_meetings/Prague/2019_10_04.LINE1")

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
library(karyoploteR)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# clean repeatMasker path
rmsk_path <- file.path(genome_dir, "rmsk.mm10.20180919.joined_rmsk_id.fa.out.gz")

# reduced exons path
exons_path <- file.path(genome_dir, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.reducedExons.RDS")

# gene info path
genes_info_path <- file.path(path = genome_dir, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.geneInfo.csv")

# selected LINE1 elements
line1_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/filter_LINE1/LINE1.full_length_and_ORF.csv"

######################################################## READ DATA
# read raw repeatMasker
rmsk_df <- readr::read_delim(file = rmsk_path, delim = "\t")

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read selected LINE1 elements
line1_tb <- readr::read_csv(line1_path)

######################################################## MAIN CODE
### get ranges of protein coding genes
# get list of protein coding
protein_coding <- 
  genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# clean exons
exons_gr %<>% 
  .[names(.) %in% protein_coding] %>% 
  range(.) %>% 
  unlist(.)

### clean repeatMasker
rmsk_tidy <- 
  rmsk_df %>%
  dplyr::mutate(strand = ifelse(strand %in% c("+", "-"), strand, "+")) %>% 
  dplyr::filter(repFamily == "L1") %>%
  # dplyr::filter(!str_detect(seqnames, "GL|JH|random")) %>% 
  GRanges(.)


### visualize on chromosomes
# prepare data
subset <- list(1:7, 8:14, 15:21)
chrs <- c(str_c("chr", 1:19), "chrX", "chrY")

# create table for karotype
rmsk_karyo <-
  rmsk_tidy %>%
  as_tibble(.) %>%
  dplyr::select(chr = seqnames, start, end, strand) %>%
  dplyr::mutate(y = runif(nrow(.), -0.75, 0.75)) %>%
  dplyr::group_by(chr) %>%
  # dplyr::sample_n(2000) %>%
  split(., .$strand)

# get table with selected LINE1 elements
rmsk_karyo <- 
  line1_tb %>% 
  dplyr::select(chr = seqnames, start, end, strand) %>%
  dplyr::mutate(y = runif(nrow(.), -0.75, 0.75)) %>%
  dplyr::group_by(chr) %>%
  # dplyr::sample_n(2000) %>%
  split(., .$strand)
  
### plot dots
for(n in 1:length(subset)){

  ### open file
  png(file.path(outpath, str_c("mm10.LINE1_genome_insertions.dots.selected.", n, ".png")), width = 1000, height = 1000, units = "px")

  # create basic plot
  kp <- plotKaryotype(genome = "mm10", plot.type = 2, chromosomes = chrs[subset[[n]]], cex = 3)

  # add plus strand
  kpDataBackground(kp, r0 = 0, r1 = 0.75, color = "#e6e6e6", data.panel = 1)
  kpPoints(kp, chr = rmsk_karyo[[1]]$chr, x = rmsk_karyo[[1]]$start, y = rmsk_karyo[[1]]$y,
           ymin = -1, ymax = 1, col = "red", pch = ".", cex = 5,
           r0 = 0, r1 = 0.75, data.panel = 1)

  # add minus strand
  kpDataBackground(kp, r0 = 0, r1 = 0.75, color = "#e6e6e6", data.panel = 2)
  kpPoints(kp, chr = rmsk_karyo[[2]]$chr, x = rmsk_karyo[[2]]$start, y = rmsk_karyo[[2]]$y,
           ymin = -1, ymax = 1, col = "red", pch = ".", cex = 5,
           r0 = 0, r1 = 0.75, data.panel = 2)

  ### close file
  dev.off()
}




### plot density
for(n in 1:length(subset)){
  
  png(file.path(outpath, str_c("mm10.LINE1_genome_insertions.density.", n, ".png")), width = 1000, height = 1000, units = "px",)
  kp <- plotKaryotype(genome = "mm10", plot.type = 2, chromosomes = chrs[subset[[n]]], cex = 3)
  
  kpPlotDensity(kp, rmsk_tidy[strand(rmsk_tidy) == "+"], 
                window.size = 10e5,
                r0 = 0, r1 = 0.75, data.panel = 1)
  kpPlotDensity(kp, exons_gr[strand(exons_gr) == "+"], 
                window.size = 10e5,
                col="#ddaacc", 
                r0 = 0, r1 = 0.75, data.panel = 1)
  
  kpPlotDensity(kp, rmsk_tidy[strand(rmsk_tidy) == "-"], 
                window.size = 10e5,
                r0 = 0, r1 = 0.75, data.panel = 2)
  kpPlotDensity(kp, exons_gr[strand(exons_gr) == "-"], 
                window.size = 10e5,
                col="#ddaacc", 
                r0 = 0, r1 = 0.75, data.panel = 2)
  dev.off()
  
}