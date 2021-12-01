#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other) in chunks, sum in loop
### DATE: 30. 06. 2018.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/test/class_bam")

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
  
  # initialize vector to hold number of reads in each category
  count_sums <- c("rRNA" = 0, "repeats" = 0, "exon" = 0, "other" = 0)
  
  ### connect to bam
  # open connection to bam file in chunks
  bamfile <- BamFile(bam_path, yieldSize = yield)
  open(bamfile)
  
  # set last read name placeholder
  last_reads <- NULL
  
  # load chunks of reads from bam file and classify each alignment
  while(length(chunk <- readGAlignmentsList(bamfile, param = ScanBamParam(what = "qname", flag = scanBamFlag(isFirstMateRead = isFirstInPair))))) {
    
    # unlist, set names of read
    chunk <- unlist(chunk)
    names(chunk) <- mcols(chunk)$qname
    
    # add last reads from previous chunk to begining of new chunk
    if(!is.null(last_reads)){
      chunk <- c(last_reads, chunk)
    }
    
    # get last reads with unique name in chunk
    last_reads <- chunk[names(chunk) == names(chunk)[length(chunk)]]
    
    # remove last reads from chunk
    chunk <- chunk[-which(names(chunk) %in% names(last_reads))]
    
    # transform to grglist (which gets ranges of only alignment part), unlist to GRanges
    chunk <- 
      GenomicRanges::grglist(chunk) %>% 
      unlist(.)
    
    # loop over bam file and set class to each alignment
    for(overlap_subject in names(ordered_filters)){
      
      # overlap reads
      overlap_reads <-
        IRanges::subsetByOverlaps(x = chunk, ranges = ordered_filters[[overlap_subject]], type = "within", ignore.strand = T) %>%
        names(.) %>%
        unique(.)
      
      # add reads to intialized vector if there is overlap with class
      if(length(overlap_reads) > 0){
        
        # add to count
        count_sums[overlap_subject] <- count_sums[overlap_subject] + length(overlap_reads)
        
        # filter reads out
        chunk <- chunk[!(names(chunk) %in% overlap_reads)]
        
      }
      
    }
    
    # set alignements which are not overlaping any category as "other"
    if(length(chunk) > 0){
      
      # add to count
      count_sums["other"] <- count_sums["other"] + length(unique(names(chunk)))
      
    }
    
  }
  
  # close connection to .bam
  close(bamfile)
  
  ### classify last read
  # transform to grglist (which gets ranges of only alignment part)
  last_reads <- GenomicRanges::grglist(last_reads)
  
  # unlist to GRanges
  last_reads <-
    last_reads %>%
    unlist(.)
  
  # overlap with all features, get which overlap is hierarchically highest 
  hierarchy_number <- findOverlaps(GRangesList(rmsk_rRNA, rmsk_gr, exons_gr), last_reads, ignore.strand = T) %>% queryHits %>% min
  
  # add to final count
  if(!is.infinite(hierarchy_number)){
    
    # if there is hit, add to highest in hierarchy 
    count_sums[hierarchy_number] <- count_sums[hierarchy_number] + 1
    
  }else{
    
    # if there is no hit, add to other category
    count_sums["other"] <- count_sums["other"] + 1
    
  }
 
  
  # return vector with read names
  return(count_sums)
  
}

######################################################## PATH VARIABLES
# set outpath, get paths for gtf and bam files
outpath <- getwd()

# get paths of genome reference files
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# get path of original bam file
bam_original_path <- "/common/WORK/fhorvat/test/class_bam/s_Lnc1Het_r2.SE.genome.Aligned.sortedByName.out.bam"

# get paths of reduced exons and repeatMasker
exons_path <- list.files(path = genome_dir, pattern = "ensembl.91.*reducedExons.RDS$", full.names = T)
rmsk_path <- list.files(path = genome_dir, pattern = "rmsk\\..*clean.*.gz$", full.names = T, recursive = T)

# get name of the samples
bam_name <-
  basename(bam_original_path) %>%
  stringr::str_replace_all(., ".genome.Aligned.sortedByName.out.bam", "")

# get path of merged bam (if there is one)
bam_merged_path <-
  list.files(path = outpath, pattern = ".genome.merged.*.bam$", full.names = T) %>%
  .[str_detect(string = ., pattern = bam_name)]

# determine whether .bam comes from pair-end sequencing
pair_end_seq <- str_detect(string = bam_name, pattern = "PE")

######################################################## READ DATA
# read gtf
exons_gr <- 
  readRDS(file = exons_path) %>% 
  unlist(.)

# read repeatMasker
rmsk_gr <-
  read_delim(file = rmsk_path, delim = "\t", col_types = cols(start = col_double(), end = col_double())) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

# filter repeatMasker - rRNA
rmsk_rRNA <- rmsk_gr[rmsk_gr$repClass == "rRNA"]

# set filter order
ordered_filters <-
  list(rmsk_rRNA, rmsk_gr, exons_gr) %>%
  magrittr::set_names(c("rRNA", "repeats", "exon"))

######################################################## MAIN CODE
# classify reads in bam file
reads_class_counts <- classReads(bam_path = bam_original_path, yield = 1000000, isFirstInPair = ifelse(pair_end_seq, T, NA))

# read second reads in pair if experiment was paired end
if(pair_end_seq){
  
  # class reads
  reads_class_counts_second <- classReads(bam_path = bam_original_path, yield = 1000000, isFirstInPair = F)

}

# classify reads in merged bam file if one exists
if(length(bam_merged_path) > 0){
  
  # class reads
  reads_class_counts_merged <-classReads(bam_path = bam_merged_path, yield = 1000000, isFirstInPair = NA) 
  
}

# add number of reads in each class to table, save
reads_class_sum_final <-
  tibble(read_group = names(reads_class_counts), count = reads_class_counts) %>%
  tidyr::spread(., read_group, count) %>%
  dplyr::mutate(sample_id = bam_name,
                total = rRNA + repeats + exon + other,
                genome.mapped_minus_rDNA = total - rRNA) %>%
  dplyr::select(sample_id, rRNA, repeats, exon, other, total, genome.mapped_minus_rDNA) %T>%
  readr::write_delim(., path = file.path(outpath, str_c(bam_name, ".read_stats.fast.txt")), delim = "\t")

