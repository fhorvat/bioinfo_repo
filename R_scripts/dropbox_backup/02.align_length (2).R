### INFO: gets histogram of length of read alignments in one bam file
### DATE: Wed Jul 18 14:32:18 2018
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Analysis/alignment_length")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tidyr)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS
widthReads <- function(bam_path, yield = 1000000, isFirstInPair = NA){
  
  # get number of alignments in bam file
  read_number <-
    Rsamtools::countBam(file = bam_path, param = ScanBamParam(flag = scanBamFlag(isFirstMateRead = isFirstInPair))) %$%
    records
  
  # initialize read class vector
  read_length <- rep(0, read_number)
  names(read_length) <- rep("nameplaceholder", read_number)
  
  # open connection to bam file in chunks
  bamfile <- BamFile(bam_path, yieldSize = yield)
  open(bamfile)
  
  # load chunks of reads from bam file and classify each alignment
  while(length(chunk <- readGAlignmentsList(bamfile, param = ScanBamParam(what = "qname", flag = scanBamFlag(isFirstMateRead = isFirstInPair))))) {
    
    # unlist, set names of reads, transform to grglist (which gets ranges of only alignment part)
    chunk <- unlist(chunk)
    names(chunk) <- mcols(chunk)$qname
    chunk <- 
      GenomicRanges::grglist(chunk) %>%
      unlist(.)
    
    # unlist to GRanges, get widths of alignments
    chunk_width <-
      chunk %>% 
      width(.)
    
    # find first element which is not filled
    last_element <- 
      which(read_length == 0) %>% 
      min(.)
    
    # set range in intialized vector
    if((last_element - 1 + yield <= read_number)){
      length_range <- last_element:(last_element - 1 + yield)
    } else{
      length_range <- last_element:read_number
    }
    
    # create named vector (read names and class as name)
    read_length[length_range] <- chunk_width
    names(read_length)[length_range] <- names(chunk)
    
  }
  
  # close connection to .bam
  close(bamfile)
  
  # return vector with read names
  return(read_length)
  
}

######################################################## PATH VARIABLES
# outpath
outpath <- getwd()

# inpath
inpath <- getwd()

# get bam file path, name, experiment
bam_path <- "%BAM_PATH"
bam_name <- basename(bam_path) %>% str_remove(., ".SE.genome.Aligned.sortedByCoord.out.bam")
experiment_name <- "Eliska_mESC_MosIR"

######################################################## READ DATA

######################################################## MAIN CODE
# get widths of all alignments in bam file
read_widths <- 
  widthReads(bam_path = bam_path) %>% 
  .[. != 0]

# get unique names
reads_widths_unique <- unique(names(read_widths))

# create sum table
reads_widths_sum <- tibble(read_width = 1:100, count = 0) 

# set loop step
loop_step <- 100000

# set loop index
loop_index <- seq(0, length(reads_widths_unique), by = loop_step)
if(loop_index[length(loop_index)] < length(reads_widths_unique)){
  loop_index <- c(loop_index, length(reads_widths_unique))
}

for(n in 1:(length(loop_index) - 1)){
  
  # subset vector
  reads_widths_subset <- read_widths[names(read_widths) %in% reads_widths_unique[(loop_index[n] + 1):loop_index[n + 1]]]
  
  # create tibble from vector
  reads_widths_df <- 
    tibble(read_id = names(reads_widths_subset), read_width = reads_widths_subset) %>% 
    dplyr::group_by(read_id) %>% 
    dplyr::summarize(read_width = max(read_width)) %>% 
    dplyr::group_by(read_width) %>% 
    dplyr::summarise(count = n())
  
  # join with sum vector
  reads_widths_sum %<>%  
    dplyr::left_join(., reads_widths_df, by = "read_width") %>% 
    dplyr::mutate(count.y = replace(count.y, is.na(count.y), 0), 
                  count = count.x + count.y) %>% 
    dplyr::select(read_width, count)
  
}

reads_widths_sum %<>% 
  dplyr::filter(count > 0) %>% 
  magrittr::set_colnames(x = ., value = c("read_length", "size")) %>% 
  dplyr::mutate(sample = bam_name,
                experiment = experiment_name) %T>% 
  readr::write_csv(x = ., path = file.path(outpath, str_c("align_length.hist", experiment_name, bam_name, "csv", sep = ".")))


