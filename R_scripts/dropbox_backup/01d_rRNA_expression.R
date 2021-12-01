#!/home/students/fhorvat/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other)
### DATE: 03. 09. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
# setwd("%OUTPATH")
# setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mesAur1/read_stats")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tibble)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)
library(BiocParallel)

######################################################## PATH VARIABLES
repeatmasker_path <- "/common/WORK/fhorvat/reference/golden_hamster/mesAur1/mesAur1.RepeatMasker.Refseq.txt"
outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mesAur1/read_stats"
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mesAur1"
bam_list <- list.files(path = bam_path, pattern = "*.bam$", full.names = T)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "GffToGRanges.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
# repeatMasker
rmsk_rRNA <-
  readr::read_delim(file = repeatmasker_path, delim = "\t", col_names = F) %>%
  magrittr::set_colnames(c("seqnames", "start", "end", "strand", "repName", "repClass_repFamily")) %>% 
  tidyr::separate(repClass_repFamily, into = c("repClass", "repFamily"), sep = "/") %>% 
  dplyr::filter(repClass == "rRNA") %>% 
  dplyr::mutate(strand = replace(strand, strand == "C", "*")) %>% 
  dplyr::mutate(full_pos = str_c(seqnames, ":", start, "-", end, "|", strand, "|", repClass, "|", repName)) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

######################################################## MAIN CODE
### get count of reads over rRNA
bamfiles <- Rsamtools::BamFileList(bam_list, yieldSize = 2000000)
BiocParallel::register(BiocParallel::MulticoreParam())
se_rptmsk <- GenomicAlignments::summarizeOverlaps(features = rmsk_rRNA,
                                                  reads = bamfiles,
                                                  mode = "IntersectionStrict",
                                                  singleEnd = FALSE,
                                                  ignore.strand = TRUE)

# data frame
se_df <- 
  assay(se_rptmsk) %>% 
  as.data.frame(.) %>% 
  as.tibble(.) %>% 
  magrittr::set_colnames(str_replace(colnames(.), ".Aligned.sortedByCoord.out.bam", "")) %>% 
  dplyr::mutate(full_pos = rmsk_rRNA$full_pos) %>% 
  dplyr::select(full_pos, everything()) %>% 
  dplyr::arrange(desc(s_GV_Hamster_r2))
