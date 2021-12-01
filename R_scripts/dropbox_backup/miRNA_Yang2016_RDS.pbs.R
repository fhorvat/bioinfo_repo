#!/home/students/fhorvat/R/bin/Rscript
### INFO: reads .bam files, gets position of smallRNA clusters 21-23 bp long, summarizeOverlaps over those positions, saves assay data.frame as .RDS 
### DATE: 18. 07. 2017.  
### AUTHOR: Filip Horvat
### qsub -V -q MASTER -l select=1:mem=100gb -N Yang -j oe -o /common/WORK/fhorvat/Projekti/Svoboda/Analyses/smallRNA//common/WORK/fhorvat/Projekti/Svoboda/Analyses/smallRNA/Yang_2016_PRJNA257532 miRNA_Yang2016_RDS.pbs.R

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/smallRNA/Yang_2016_PRJNA257532")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tibble)

library(Rsamtools)
library(Biostrings)
library(GenomicRanges)
library(GenomicAlignments)

######################################################## PATH VARIABLES
outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/smallRNA/Yang_2016_PRJNA257532"
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/Yang_2016_development_smallRNA_PRJNA257532/Data/Mapped/STAR_mm10"
bam_list <- list.files(path = bam_path, pattern = "s_oocyte.*bam$|s_1cell.*bam$", full.names = T)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## READ DATA

######################################################## MAIN CODE
### read coordinates of mapped reads from .bam
# create list of genomic coordinates 21-23 bp long
bam_gr_all <- 
  lapply(bam_list, function(bam_file){
    
    # read bam file as data frame
    param <- Rsamtools::ScanBamParam(what = c("rname", "strand", "pos", "seq", "cigar"))
    bam_in <- Rsamtools::scanBam(bam_file, param = param)
    
    # get aligned sequence
    seq_aligned_width <-
      GenomicAlignments::sequenceLayer(x = bam_in[[1]]$seq, cigar = bam_in[[1]]$cigar) %>%
      width(.)
    
    # convert to data.frame, filter for miRNA clusters
    bam_in[[1]]$seq <- NULL
    bam_in[[1]]$cigar <- NULL
    bam_gr <-
      dplyr::bind_cols(bam_in) %>%
      dplyr::mutate(seq_width = seq_aligned_width) %>%
      dplyr::mutate(end = (pos + seq_width - 1),
                    full_pos = str_c(rname, ":", pos, "-", end, "|", strand)) %>%
      dplyr::filter(!duplicated(full_pos)) %>%
      dplyr::select(-seq_width, seqnames = rname, start = pos, end, strand, full_pos) %>%
      makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
      GenomicRanges::reduce(.) %>%
      .[width(.) >= 21 & width(.) <= 23]
    
    return(bam_gr)
    
  }) %>% 
  Reduce(function(gr1, gr2) c(gr1, gr2), .) %>% 
  unique(.) %>% 
  GenomicRanges::reduce(.) %>% 
  .[width(.) >= 21 & width(.) <= 23]


### counts and RPMs
# summarizeOverlaps
bamfiles <- Rsamtools::BamFileList(bam_list, yieldSize = 2000000)
BiocParallel::register(BiocParallel::MulticoreParam())
se <- GenomicAlignments::summarizeOverlaps(features = bam_gr_all,
                                           reads = bamfiles,
                                           mode = "IntersectionStrict",
                                           singleEnd = TRUE,
                                           ignore.strand = FALSE)

# get assay data.frame, gather to long data.frame
assay_df <- 
  assay(se) %>% 
  as_tibble(.) %>% 
  magrittr::set_colnames(., str_replace(colnames(.), "_Aligned.sortedByCoord.out.bam", "")) %>% 
  dplyr::mutate(full_pos = str_c(seqnames(bam_gr_all), ":", start(bam_gr_all), "-", end(bam_gr_all), "|", strand(bam_gr_all))) %>% 
  tidyr::gather(key = sample, value = count, -full_pos)

# save as .RDS
saveRDS(object = assay_df, file = file.path(outpath, "smallRNA_reduced_assay_df.RDS"))