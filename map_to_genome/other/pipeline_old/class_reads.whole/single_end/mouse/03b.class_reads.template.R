#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other)
### DATE: 04. 10. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("%OUT_PATH")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# counts reads in each category (rDNA -> repeat -> exon -> other) with subsequent filtering
classReads <- function(granges_bam, pair_end = pair_end_seq){
  
  # set filter order
  ordered_filters <-
    list(rmsk_rRNA, rmsk_gr, refseq_exons) %>%
    magrittr::set_names(c("rRNA", "repeats", "exon"))
  reads_sum <- tibble(sample_id = bam_name)
  
  # loop over bam file
  for(overlap_subject in names(ordered_filters)){
    
    # overlap reads
    overlap_reads <-
      IRanges::subsetByOverlaps(x = granges_bam, ranges = ordered_filters[[overlap_subject]], type = "within", ignore.strand = T) %>%
      names(.) %>%
      unique(.)
    
    # add columnt to sum table
    reads_sum %<>%
      mutate(overlap_subject = length(overlap_reads)) %>%
      data.table::setnames(., old = ncol(.), new = overlap_subject)
    
    # filter reads out
    granges_bam <- granges_bam[!(names(granges_bam) %in% overlap_reads)]
    
  }
  
  # set everything else to "other" and sum total
  reads_sum %<>%
    dplyr::mutate(other = length(unique(names(granges_bam)))) %>%
    dplyr::mutate(total = rowSums(.[, -1]), 
                  genome.mapped_minus_rDNA = total - rRNA, 
                  genome.mapped_minus_rDNA = ifelse(pair_end, round(genome.mapped_minus_rDNA / 2), genome.mapped_minus_rDNA))
  
  return(reads_sum)
  
}

######################################################## PATH VARIABLES
# set outpath, get paths for gtf and bam files
outpath <- getwd()

# get paths of refSeq and repeatMasker
genome_dir <- "/common/WORK/fhorvat/reference/mouse/mm10"
refseq_path <- list.files(path = genome_dir, pattern = "refGene.*gtf.gz$", full.names = T, recursive = T)
rmsk_path <- 
  list.files(path = genome_dir, pattern = "rmsk.*clean.*.gz$", full.names = T, recursive = T) %>% 
  ifelse(test = any(str_detect(., "VIZ")), yes = .[str_detect(., "VIZ")], no = .)

# get path of original bam file
bam_original_path <- "%BAM_PATH"

# get name of the samples
bam_name <-
  basename(bam_original_path) %>%
  stringr::str_replace_all(., ".genome.Aligned.sortedByCoord.out.bam", "")

# get path of merged bam (if there is one)
bam_merged_path <- 
  list.files(path = outpath, pattern = ".genome.merged.", full.names = T) %>% 
  .[str_detect(string = ., pattern = bam_name)]

# determine whether .bam comes from pair-end sequencing
pair_end_seq <- str_detect(string = bam_name, pattern = "PE")

######################################################## READ DATA
# refSeq gtf
refseq_exons <-
  read_delim(file = refseq_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c")) %>%
  gtfToGRanges(., "exon")

# repeatMasker
rmsk_gr <-
  read_delim(file = rmsk_path, delim = "\t", col_types = cols(start = col_double(), end = col_double())) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

# repeatMasker - rRNA
rmsk_rRNA <- rmsk_gr[rmsk_gr$repClass == "rRNA"]

######################################################## MAIN CODE
### count reads over features hierarchically - complete match: rDNA -> repeat -> exon -> other
# read pair-end or single end bam file, if pair-end set unique names for each mate in fragment 
if(pair_end_seq){
  gbam_original <- bamToGRangesList(bam_original_path)
  names(gbam_original) <- str_c(names(gbam_original), 1:2)
  gbam_original <- 
    unlist(gbam_original) %>% 
    classReads(granges_bam = .)
}else{
  gbam_original <-
    bamToGRangesList(bam_original_path) %>%
    unlist(.) %>% 
    classReads(granges_bam = .)
}

# read merged bam if there is one
if(length(bam_merged_path) > 0){
  gbam_merged <-
    bamToGRangesList(bam_merged_path) %>%
    unlist(.) %>% 
    classReads(granges_bam = .)
}else{
  gbam_merged <- NULL
}

# merge both, write
rbind(gbam_original, gbam_merged) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::summarise_all(.funs = sum) %>% 
  readr::write_delim(., path = file.path(outpath, str_c(bam_name, ".read_stats.txt")), delim = "\t")
