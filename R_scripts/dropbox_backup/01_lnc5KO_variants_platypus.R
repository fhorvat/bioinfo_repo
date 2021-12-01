#!/home/students/fhorvat/R/bin/Rscript
### INFO: read VCF from lnc5 KOs and WT, find exons with mutations in lnc5 KO which are not mutated in WT
### DATE: 22. 5. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/variant_calling")

######################################################## LIBRARIES
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(magrittr)
library(readr)
library(stringr)
library(tibble)
library(data.table)

library(GenomicFeatures)
library(VariantAnnotation)

######################################################## PATH VARIABLES
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/variant_calling/platypus/vcf_output"

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/variant_calling/platypus/results"

lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"

gtf_path_ENSEMBL <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Ensembl_GRCm38.86.20161128.gtf.gz"

######################################################## SOURCE FILES
source(file.path(lib_path, "GffToGRanges.R"))

######################################################## FUNCTIONS
# .vcf to ranges with filtering
readFilterVCF <- function(vcf_path, min_coverage, filter_pass = F){
  
  # read .vcf file
  vcf <- readVcf(file = vcf_path)
  
  # get ranges, add metadata
  ranges_vcf <- rowRanges(vcf)
  mcols(ranges_vcf)$NV <- geno(vcf)$NV
  ranges_vcf <- ranges_vcf[sapply(mcols(ranges_vcf)$NV, length) == 1]
  mcols(ranges_vcf)$NV <- unlist(mcols(ranges_vcf)$NV)
  
  # filter ranges
  if(filter_pass){
    ranges_vcf <- ranges_vcf[mcols(ranges_vcf)$FILTER == "PASS"]
  }
  ranges_vcf <- ranges_vcf[mcols(ranges_vcf)$NV >= min_coverage]
  
  return(ranges_vcf)
  
}

# finds overlaps between query and subject, filters query overlaping the subject
filterByNames <- function(query, subject){
  query_filtered <- query[!(names(query) %in% names(subject))]
  return(query_filtered)
}

# finds intersect between names of query and subject, returns intersect
intersectByNames <- function(query, subject){
  query_filtered <- query[names(query) %in% names(subject)]
  return(query_filtered)
}

######################################################## READ DATA
### read .gtf
# ENSEMBL gtf
gtf <-
  read_delim(file = gtf_path_ENSEMBL, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c")) %>% 
  GffToGRanges(., filter = "exon") 

### read .vcf
# get list of .vcf files
vcf_list <- list.files(path = inpath, pattern = "*.vcf.gz$", full.names = T)

# lnc5 KO mutated indels
ranges_lnc5 <- lapply(vcf_list[str_detect(vcf_list, "Lnc5")], readFilterVCF, min_coverage = 3, filter_pass = F)

# WT mutated indels
ranges_WT <- lapply(vcf_list[str_detect(vcf_list, "WT")], readFilterVCF, min_coverage = 3, filter_pass = F)

######################################################## MAIN CODE
### keep mutations present in all lnc5 KO replicates 
ranges_lnc5_intersect <- Reduce(f = intersectByNames, x = ranges_lnc5)

### filter out WT mutations from lnc5 KO mutations
invisible(lapply(c("ALL", "ANY"), function(filter){
  
  if(filter == "ALL"){
    # filter those which appear in ALL WT replicates
    ranges_WT_intersect <- Reduce(f = intersectByNames, x = ranges_WT)
    ranges_lnc5_filtered <- filterByNames(query = ranges_lnc5_intersect, subject = ranges_WT_intersect)
  }else{
    if(filter == "ANY"){
      # filter those which appear in ANY WT replicates
      ranges_lnc5_filtered <- Reduce(f = filterByNames, x = ranges_WT, init = ranges_lnc5_intersect)
    }
  }
  
  # take only mismatches in exons
  ranges_lnc5_filtered <- unique(ranges_lnc5_filtered[queryHits(findOverlaps(ranges_lnc5_filtered, gtf))])
  
  # mismatch loci output
  ranges_lnc5_filtered_df <- 
    ranges_lnc5_filtered %>% 
    as.data.frame() %>% 
    dplyr::select(-c(REF, ALT, QUAL, width, strand, paramRangeID, FILTER), mis_coverage = NV) %>% 
    tibble::rownames_to_column(var = "mismatch") %>% 
    dplyr::mutate(mismatch = str_replace(mismatch, ".*_", "")) %>%
    tidyr::separate(mismatch, c("ref", "mismatch"), "/") %>%
    dplyr::select(3:5, 1:2, 6) %T>% 
    readr::write_csv(x = ., path = file.path(outpath, str_c("lnc5KO_mismatches_", filter, "_platypus.csv")))
  
  # mismatch exon output
  mismatch_exon_overlaps <- findOverlaps(gtf, ranges_lnc5_filtered)
  mutated_exons <- 
    gtf[queryHits(mismatch_exon_overlaps)] %>% 
    as.data.frame() %>% 
    dplyr::mutate(mismatch = names(ranges_lnc5_filtered[subjectHits(mismatch_exon_overlaps)]), 
                  NV = mcols(ranges_lnc5_filtered[subjectHits(mismatch_exon_overlaps)])$NV, 
                  fullName = str_c(strand, mismatch)) %>% 
    dplyr::filter(!duplicated(fullName)) %>% 
    dplyr::select(-c(width, frame, feature.type, fullName, transcript_id, exon_id), mis_coverage = NV) %>% 
    dplyr::mutate(mismatch = str_replace(mismatch, "chr[1-9,MT].*:", "")) %>%
    tidyr::separate(mismatch, c("mismatch_pos", "ref", "mismatch"), "_|/") %T>%
    readr::write_csv(x = ., path = file.path(outpath, str_c("lnc5KO_mismatches_exons_", filter, "_platypus.csv")))
  
}))