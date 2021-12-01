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

library(GenomicFeatures)
library(VariantAnnotation)

######################################################## PATH VARIABLES
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/variant_calling/GATK/vcf_output"

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/variant_calling/GATK/results"

lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"

gtf_path_ENSEMBL <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Ensembl_GRCm38.86.20161128.gtf.gz"

######################################################## SOURCE FILES
source(file.path(lib_path, "GffToGRanges.R"))

######################################################## FUNCTIONS
# .vcf to ranges with filtering
readFilterVCF <- function(vcf_path, min_coverage){
  
  # read .vcf file
  vcf <- readVcf(file = vcf_path)
  
  # get ranges 
  ranges_vcf <- 
    do.call(cbind, geno(vcf)@listData) %>% 
    as.data.frame() %>% 
    set_colnames(c("genotype", "allele_coverage", "total_coverage", "quality", "phred")) %>% 
    tibble::rownames_to_column(., var = "mismatch_name") %>% 
    as_tibble() %>% 
    cbind(., as_tibble(do.call(rbind, .$allele_coverage))) %>% 
    data.table::setnames(., c("V1", "V2", "V3"), c("ref_coverage", "mismatch_coverage", "mismatch2_coverage")) %>% 
    tidyr::separate(mismatch_name, c("seqnames", "mismatch"), ":", remove = F) %>% 
    tidyr::separate(mismatch, c("start", "ref", "mismatch"), "_|/") %>% 
    dplyr::mutate(end = start) %>% 
    dplyr::select(seqnames, start, end, ref, mismatch, ref_coverage, mismatch_coverage, quality, mismatch_name) %>% 
    dplyr::filter(mismatch_coverage >= min_coverage) %>%
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T) %>% 
    set_names(., mcols(.)$mismatch_name)
    
  # remove mismatch_name from mcols
  mcols(ranges_vcf)$mismatch_name <- NULL
  
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
# ENSEMBL gtf
gtf <-
  read_delim(file = gtf_path_ENSEMBL, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c")) %>% 
  GffToGRanges(., filter = "exon") 

# get list of .vcf files
vcf_list <- list.files(path = inpath, pattern = "*.vcf.gz$", full.names = T)

# lnc5 KO mutated indels
ranges_lnc5 <- lapply(vcf_list[str_detect(vcf_list, "Lnc5")], readFilterVCF, min_coverage = 3)

# WT mutated indels
ranges_WT <- lapply(vcf_list[str_detect(vcf_list, "WT")], readFilterVCF, min_coverage = 3)

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
  
  ### mismatch loci output
  ranges_lnc5_filtered_df <- 
    ranges_lnc5_filtered %>% 
    as.data.frame() %>% 
    set_rownames(., NULL) %>% 
    dplyr::select(-c(width, strand, quality)) %>% 
    readr::write_csv(x = ., path = file.path(outpath, str_c("lnc5KO_mismatches_", filter, "_GATK.csv")))
  
  ### mismatch exons output
  mismatch_exon_overlaps <- findOverlaps(gtf, ranges_lnc5_filtered)
  mutated_exons_df <-
    gtf[queryHits(mismatch_exon_overlaps)] %>%
    as.data.frame() %>%
    dplyr::mutate(mismatch = names(ranges_lnc5_filtered[subjectHits(mismatch_exon_overlaps)]),
                  ref_coverage = mcols(ranges_lnc5_filtered[subjectHits(mismatch_exon_overlaps)])$ref_coverage,
                  mismatch_coverage = mcols(ranges_lnc5_filtered[subjectHits(mismatch_exon_overlaps)])$mismatch_coverage,
                  fullName = str_c(strand, mismatch)) %>%
    dplyr::filter(!duplicated(fullName)) %>%
    dplyr::select(-c(width, frame, feature.type, fullName, gene_biotype, exon_id)) %>%
    dplyr::mutate(mismatch = str_replace(mismatch, "chr[1-9,MT].*:", "")) %>%
    tidyr::separate(mismatch, c("mismatch_pos", "ref", "mismatch"), "_|/") %T>%
    readr::write_csv(x = ., path = file.path(outpath, str_c("lnc5KO_mismatches_exons_", filter, "_GATK.csv")))
  
}))
