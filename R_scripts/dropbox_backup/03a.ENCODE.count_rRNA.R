#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other) in chunks
### DATE: 21. 06. 2018.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/ENCODE_2014_Nature_GSE49417/Data/Mapped/STAR_mm10")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(data.table)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
## counts reads in each category (rDNA -> repeat -> exon -> other) with subsequent filtering
classReads <- function(bam_path, yield = 1000000, isFirstInPair = NA){
  
  # get number of alignments in bam file
  read_number <- 
    Rsamtools::countBam(file = bam_path, param = ScanBamParam(flag = scanBamFlag(isFirstMateRead = isFirstInPair, isSecondaryAlignment = F))) %$% 
    records
  
  # initialize read class vector
  read_class <- NULL

  # open connection to bam file in chunks
  bamfile <- BamFile(bam_path, yieldSize = yield)
  open(bamfile)
  
  # load chunks of reads from bam file and classify each alignment
  while(length(chunk <- readGAlignmentsList(bamfile, param = ScanBamParam(what = "qname", flag = scanBamFlag(isFirstMateRead = isFirstInPair))))) {
    
    # unlist, set names of reads, transform to grglist (which gets ranges of only alignment part)
    chunk <- unlist(chunk)
    names(chunk) <- mcols(chunk)$qname
    chunk <- GenomicRanges::grglist(chunk)
    
    # unlist to GRanges
    chunk <- 
      chunk %>% 
      unlist(.) 
    
    # overlap reads with rRNA
    overlap_reads <-
      IRanges::subsetByOverlaps(x = chunk, ranges = rmsk_rRNA, type = "within", ignore.strand = T) %>%
      names(.) %>%
      unique(.)
    
    # add reads to intialized vector if there is overlap with class
    if(length(overlap_reads) > 0){
      
      # create named vector (read names and class as name)
      read_class <- c(read_class, overlap_reads)
      
    }
    
  }
  
  # close connection to .bam
  close(bamfile)
  
  # create tibble with number of reads in bam file classifed as rRNA, set the rest to other
  read_class_df <- 
    tibble(rRNA = length(unique(read_class)), 
           total = read_number, 
           other = total - rRNA)
  
  # return vector with read names
  return(read_class_df)
  
}

######################################################## PATH VARIABLES
# set outpath, get paths for gtf and bam files
outpath <- getwd()

# get paths of genome reference files
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# get path of original bam file
bam_original_path <- "%BAM_PATH"

# get paths of reduced exons and repeatMasker
rmsk_path <- list.files(path = genome_dir, pattern = "rmsk\\..*clean.*.gz$", full.names = T, recursive = T)

# get name of the samples
bam_name <-
  basename(bam_original_path) %>%
  stringr::str_replace_all(., ".genome.Aligned.sortedByCoord.out.bam", "")

# get path of merged bam (if there is one)
bam_merged_path <-
  list.files(path = outpath, pattern = ".genome.merged.*.bam$", full.names = T) %>%
  .[str_detect(string = ., pattern = bam_name)]

# determine whether .bam comes from pair-end sequencing
pair_end_seq <- str_detect(string = bam_name, pattern = "PE")

######################################################## READ DATA
# read repeatMasker
rmsk_rRNA <-
  read_delim(file = rmsk_path, delim = "\t", col_types = cols(start = col_double(), end = col_double())) %>%
  dplyr::filter(repClass == "rRNA") %>% 
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

######################################################## MAIN CODE
### classify reads in bam file
# a) first read in pair
peakRAM(reads_class_first <- classReads(bam_path = bam_original_path, yield = 1000000, isFirstInPair = T))

# b) second read in pair
reads_class_second <- classReads(bam_path = bam_original_path, yield = 1000000, isFirstInPair = F) 

# c) merged bam file if one exists
if(length(bam_merged_path) > 0){
  reads_class_vector_merged <- classReads(bam_path = bam_merged_path, yield = 1000000, isFirstInPair = NA)
}
  
# sum number of reads in each class
reads_class_sum <- 
  dplyr::bind_rows(reads_class_first, 
                   reads_class_second, 
                   reads_class_vector_merged) %>% 
  dplyr::summarise_all(funs(sum(.))) %>% 
  dplyr::mutate(sample_id = bam_name, 
                repeats = 0,
                exon = 0, 
                genome.mapped_minus_rDNA = total - rRNA) %>% 
  dplyr::select(sample_id, rRNA, repeats, exon, other, total, genome.mapped_minus_rDNA) %T>%
  readr::write_delim(., path = file.path(outpath, str_c(bam_name, ".read_stats.txt")), delim = "\t")

