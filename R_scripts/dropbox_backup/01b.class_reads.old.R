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
  
  # set filter order
  ordered_filters <-
    list(rmsk_rRNA, rmsk_gr, exons_gr) %>%
    magrittr::set_names(c("rRNA", "repeats", "exon"))
  
  # get number of alignments in bam file
  read_number <-
    Rsamtools::countBam(file = bam_path, param = ScanBamParam(flag = scanBamFlag(isFirstMateRead = isFirstInPair))) %$%
    records
  
  # initialize read class vector
  read_class <- rep("placeholder", read_number)
  names(read_class) <- rep("nameplaceholder", read_number)
  
  # open connection to bam file in chunks
  bamfile <- BamFile(bam_path, yieldSize = yield)
  open(bamfile)
  
  # load chunks of reads from bam file and classify each alignment
  while(length(chunk <- readGAlignmentsList(bamfile, param = ScanBamParam(what = "qname", flag = scanBamFlag(isFirstMateRead = isFirstInPair))))) {
    
    # unlist, set names of reads, transform to grglist (which gets ranges of only alignment part)
    chunk <- unlist(chunk)
    names(chunk) <- mcols(chunk)$qname
    chunk <- GenomicRanges::grglist(chunk)
    
    # unlist to GRanges, classify alignments in loop
    chunk <-
      chunk %>%
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
        
        # find first element which is not filled
        last_element <-
          which(read_class == "placeholder") %>%
          min(.)
        
        # set range in intialized vector
        class_range <- last_element:(last_element - 1 + length(overlap_reads))
        
        # create named vector (read names and class as name)
        read_class[class_range] <- overlap_reads
        names(read_class)[class_range] <- overlap_subject
        
        # filter reads out
        chunk <- chunk[!(names(chunk) %in% overlap_reads)]
        
      }
      
    }
    
    # set alignements which are not overlaping any category as "other"
    if(length(chunk) > 0){
      
      # find first element which is not filled
      last_element <-
        which(read_class == "placeholder") %>%
        min(.)
      
      # set range in intialized vector
      class_range <- last_element:(last_element - 1 + length(unique(names(chunk))))
      
      # add to named vector (read names and class as name)
      read_class[class_range] <- unique(names(chunk))
      names(read_class)[class_range] <- "other"
      
    }
    
  }
  
  # close connection to .bam
  close(bamfile)
  
  # return vector with read names
  return(read_class)
  
}

######################################################## PATH VARIABLES
# set outpath, get paths for gtf and bam files
outpath <- getwd()

# get paths of genome reference files
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# get path of original bam file
bam_original_path <- "./s_Lnc1Het_r2.SE.genome.Aligned.sortedByCoord.out.bam"

# get paths of reduced exons and repeatMasker
exons_path <- list.files(path = genome_dir, pattern = "ensembl.91.*reducedExons.RDS$", full.names = T)
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
# read gtf
exons_gr <- readRDS(file = exons_path)

# read repeatMasker
rmsk_gr <-
  read_delim(file = rmsk_path, delim = "\t", col_types = cols(start = col_double(), end = col_double())) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

# filter repeatMasker - rRNA
rmsk_rRNA <- rmsk_gr[rmsk_gr$repClass == "rRNA"]

######################################################## MAIN CODE
# classify reads in bam file
reads_class_vector <-
  classReads(bam_path = bam_original_path, yield = 1000000, isFirstInPair = ifelse(pair_end_seq, T, NA)) %>%
  .[. != "placeholder"]

# read second reads in pair if experiment was paired end
if(pair_end_seq){
  
  # class reads
  reads_class_vector_second <-
    classReads(bam_path = bam_original_path, yield = 1000000, isFirstInPair = F) %>%
    .[. != "placeholder"]
  
  # set unique name
  reads_class_vector_second %<>%
    str_c(., ".1") %>%
    set_names(., names(reads_class_vector_second))
  
  # join with original reads
  reads_class_vector <- c(reads_class_vector, reads_class_vector_second)
  
}

# classify reads in merged bam file if one exists
if(length(bam_merged_path) > 0){
  
  reads_class_vector_merged <-
    classReads(bam_path = bam_merged_path, yield = 1000000, isFirstInPair = NA) %>%
    .[. != "placeholder"]
  
  # join with original reads
  reads_class_vector <- c(reads_class_vector, reads_class_vector_merged)
  
}

# remove from memory
rm(list = ls()[str_detect(ls(), "reads_class_vector_second|reads_class_vector_merged")]); gc()

### loop through unique names and sum read categories
# get unique names
reads_class_vector_unique <- unique(reads_class_vector)

# create sum table
reads_class_sum <- tibble(read_group = c("rRNA", "repeats", "exon", "other"), count = 0)

# set loop step
loop_step <- 100000

# set loop index
loop_index <- seq(0, length(reads_class_vector_unique), by = loop_step)
if(loop_index[length(loop_index)] < length(reads_class_vector_unique)){
  loop_index <- c(loop_index, length(reads_class_vector_unique))
}

# loop
for(n in 1:(length(loop_index) - 1)){
  
  # subset vector
  reads_class_vector_subset <- reads_class_vector[reads_class_vector %in% reads_class_vector_unique[(loop_index[n] + 1):loop_index[n + 1]]]
  
  # create tibble from vector
  reads_class_dt <- data.table(read_id = reads_class_vector_subset, read_class = names(reads_class_vector_subset))
  
  # classify reads over features hierarchically - complete match: rDNA -> repeat -> exon -> other
  reads_class_dt[ , read_group := ifelse("rRNA" %in% read_class, "rRNA",
                                         ifelse("repeats" %in% read_class, "repeats",
                                                ifelse("exon" %in% read_class, "exon",
                                                       ifelse("other" %in% read_class, "other", "no_class")))),
                  by = read_id]
  
  # summarize
  reads_class_sum_subset <-
    reads_class_dt %>%
    dplyr::select(read_id, read_group) %>%
    unique(.) %>%
    dplyr::group_by(read_group) %>%
    dplyr::summarise(count = n())
  
  # join with sum vector
  reads_class_sum %<>%
    dplyr::left_join(., reads_class_sum_subset, by = "read_group") %>%
    dplyr::mutate(count.y = replace(count.y, is.na(count.y), 0),
                  count = count.x + count.y) %>%
    dplyr::select(read_group, count)
  
}

# sum number of reads in each class
reads_class_sum_final <-
  reads_class_sum %>%
  tidyr::spread(., read_group, count) %>%
  dplyr::mutate(sample_id = bam_name,
                total = rRNA + repeats + exon + other,
                genome.mapped_minus_rDNA = total - rRNA) %>%
  dplyr::select(sample_id, rRNA, repeats, exon, other, total, genome.mapped_minus_rDNA) %T>%
  readr::write_delim(., path = file.path(outpath, str_c(bam_name, ".read_stats.txt")), delim = "\t")

