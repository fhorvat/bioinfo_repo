#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads in bam file sorted by read names into categories (rRNA, repeat, exon, other) in chunks
### DATE: Mon Mar 04 14:51:59 2019
### AUTHOR: Filip Horvat

# arguments from command line
args <- commandArgs(TRUE)

######################################################## WORKING DIRECTORY

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

  # initialize empty data.tabe
  read_counts_empty <- data.table(hits = 1:3)

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

    # get number of unique reads in chunk
    unique_reads_number <- length(unique(names(chunk)))

    # transform to grglist (which gets ranges of only alignment part), unlist to GRanges
    chunk <-
      GenomicRanges::grglist(chunk) %>%
      unlist(.)

    # count reads overlaping rDNA, repeats and exons
    read_counts <- findOverlaps(chunk, features_list, ignore.strand = T, type = "within")

    # create data.table with counts
    read_counts <-
      data.table(hits = subjectHits(read_counts),
                 reads = names(chunk[queryHits(read_counts)])) %>%
      .[, list(hits = min(hits)), by = "reads"] %>%
      .[order(hits), list(n = .N), by = "hits"] %>%
      .[read_counts_empty, on = "hits"] %>%
      .[is.na(n), n := 0] %$%
      n

    # calculate how many reads are not overlaping any class
    read_counts_other <- unique_reads_number - sum(read_counts)

    # add other reads to read counts
    read_counts <- c(read_counts, read_counts_other)

    # add to count sums
    count_sums <- count_sums + read_counts

  }

  # close connection to .bam
  close(bamfile)

  ### classify last read
  # transform to grglist (which gets ranges of only alignment part), unlist to GRanges
  last_reads <-
    GenomicRanges::grglist(last_reads) %>%
    unlist(.)

  # overlap with all features, get which overlap is hierarchically highest
  hierarchy_number <- findOverlaps(features_list, last_reads, ignore.strand = T) %>% queryHits %>% min

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
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# genome directory, gtf and repeatMasker
bam_original_path <- args$bam_file
gtf_path <- args$gtf_path
rmsk_path <- args$rmsk_path

# get name of the samples
bam_name <- basename(bam_original_path)

# check if bam is sorted by read name
if(str_detect(bam_name, "\\.sortedByName\\.bam$")){
  
  # remove sufix from bam name
  bam_name %<>% 
    stringr::str_remove_all(., "\\.sortedByName\\.bam")
  
}else{
  
  # kill script
  stop("Something's not right - did you sort .bam by read names? If not, please do, this script will produce wrong results if run on bam which is not sorted by read names!")
  
}

# determine whether .bam comes from pair-end sequencing
pair_end_seq <- str_detect(string = bam_name, pattern = "PE")

######################################################## READ DATA
# read .gtf
gtf_tb <- read_delim(file = gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# read repeatMasker
rmsk_gr <-
  read_delim(file = rmsk_path, delim = "\t", col_types = cols(start = col_double(), end = col_double())) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

# filter repeatMasker - rRNA
rmsk_rRNA <- rmsk_gr[rmsk_gr$repClass == "rRNA"]

######################################################## MAIN CODE
### get reduced exons coordinates, calculate length of exons
# convert GTF to GRanges, get only exons, reduce
exons_gr <-
  gtfToGRanges(gtf_tb, filter = "exon") %>%
  GenomicRanges::split(., .$gene_id) %>%
  GenomicRanges::reduce(., ignore.strand = F) %>%
  unlist(.)

# join to one list
features_list <- GRangesList("rRNA" = rmsk_rRNA, "repeats" = rmsk_gr, "exon" = exons_gr)

# classify reads in bam file
reads_class_counts <- classReads(bam_path = bam_original_path, yield = 1000000, isFirstInPair = ifelse(pair_end_seq, T, NA))

# read second reads in pair if experiment was paired end
if(pair_end_seq){
  
  # class reads
  reads_class_counts_second <- classReads(bam_path = bam_original_path, yield = 1000000, isFirstInPair = F)
  
  # sum with counts
  reads_class_counts <- reads_class_counts + reads_class_counts_second
  
  # get number of FRAGMENTS which are mapping to genome except rDNA (for normalizing)
  mapped_minus_rDNA <-
    (sum(reads_class_counts) - unname(reads_class_counts["rRNA"])) %>%
    magrittr::divide_by(., 2) %>%
    round(.)
  
  # add number of reads in each class to table, save
  reads_class_sum_final <-
    tibble(read_group = names(reads_class_counts), count = reads_class_counts) %>%
    tidyr::spread(., read_group, count) %>%
    dplyr::mutate(sample_id = bam_name,
                  total = rRNA + repeats + exon + other,
                  mapped_minus_rDNA = mapped_minus_rDNA) %>%
    dplyr::select(sample_id, rRNA, repeats, exon, other, total, mapped_minus_rDNA) %T>%
    readr::write_delim(., path = file.path(outpath, str_c(bam_name, ".read_stats.txt")), delim = "\t")
  
}else{
  
  # add number of reads in each class to table, save
  reads_class_sum_final <-
    tibble(read_group = names(reads_class_counts), count = reads_class_counts) %>%
    tidyr::spread(., read_group, count) %>%
    dplyr::mutate(sample_id = bam_name,
                  total = rRNA + repeats + exon + other,
                  mapped_minus_rDNA = total - rRNA) %>%
    dplyr::select(sample_id, rRNA, repeats, exon, other, total, mapped_minus_rDNA) %T>%
    readr::write_delim(., path = file.path(outpath, str_c(bam_name, ".read_stats.txt")), delim = "\t")
  
}



