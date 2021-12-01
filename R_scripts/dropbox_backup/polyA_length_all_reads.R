### INFO: R Script
### DATE: 22. 06. 2017. 
### AUTHOR: Filip Horvat
### PATH: 
rm(list = ls()); gc()
# options(bitmapType = 'cairo')

################################################################################### WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Analysis/polyA_length")

################################################################################### LIBRARIES
# data shaping
library(dplyr)
library(data.table)
library(stringr)
library(readr)
library(magrittr)
library(ggplot2)
library(reshape2)
library(tibble)

# genomics
library(ShortRead)

################################################################################### PATH VARIABLES
fastq_path <- "/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Data/Raw/Cleaned/02_demultiplexed/bbtools_seal"
fastq_files <- list.files(fastq_path, pattern = "[B6|DBA]_[GV|MII].*fastq$", full.names = T)

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Analysis/polyA_length"
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"

################################################################################### SOURCE FILES

################################################################################### FUNCTIONS

################################################################################### SCRIPT PARAMS

################################################################################### TABLES 

################################################################################### MAIN CODE
# loop through fastq files
all_reads_polyA <- lapply(1:length(fastq_files), function(X){
  
  # read fastq file, get sequences and read ID's
  fastq_all <- ShortRead::readFastq(dirPath = fastq_files[X])
  fastq_seq <- sapply(ShortRead::sread(fastq_all), toString)
  
  ### problem: 
  # - fuzzy matching function in R (aregexec) returns only first pattern occurance 
  # - if polyA tail length is set too low it won't find longer polyA sequences in same read downstream
  # solution:
  # - set polyA tail length to some big value (~25) and find all sequences with polyA of that length or longer
  # - sequentially lower polyA tail length by 1 and do pattern matching again with lower value
  # - repeat until you have longest polyA sequence in all reads
  
  ## initialize values
  # full fastq table with all polyA length values (set inital polyA length to NA)
  fastq_table_full <- tibble(name = 1:length(fastq_seq),
                             polyA_final_length = NA, 
                             polyA_final_pos = NA, 
                             fastq_file = str_replace_all(string = fastq_files[X], pattern = "\\/.*\\/|.fastq", replacement = ""), 
                             fastq_seq = fastq_seq)
  
  # temporary table with added sequences 
  fastq_table <- fastq_table_full
  
  # initial polyA length
  polyA_min_length <- 25
  
  # loop - sequentially lower polyA tail length by 1 
  repeat{
    
    # keep only NA polyA lengths (so pattern matching is done only on sequences without polyA length value defined before)
    fastq_table %<>% 
      dplyr::filter(is.na(polyA_final_length))
    
    # polyA fuzzy pattern match 
    polyA_match <- aregexec(pattern = str_c("A{", polyA_min_length, ",}|T{", polyA_min_length, ",}"), 
                            text = fastq_table$fastq_seq, 
                            max.distance = 0.1)
      
    # polyA position
    polyA_pos <- 
      polyA_match %>% 
      unlist(.)
    
    # get lengths of polyA tails in sequences from temporary table (if there is no match for that polyA length set value to 0)
    polyA_length <-
      regmatches(fastq_table$fastq_seq, polyA_match) %>%
      lapply(., nchar) %>% 
      lapply(., function(x) replace(x, length(x) == 0, 0)) %>% 
      unlist(.)
    
    # break loop if polyA length vector is empty
    if(length(polyA_length) == 0){
      break()
    }
    
    # replace all polyA lengths > 0
    fastq_table %<>%
      dplyr::select(name, polyA_final_length, polyA_final_pos, fastq_seq) %>% 
      tibble::add_column(polyA_length = polyA_length, 
                         polyA_pos = polyA_pos) %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(polyA_final_length = replace(polyA_final_length, polyA_length != 0, polyA_length), 
                    polyA_final_pos = replace(polyA_final_pos, polyA_pos != -1, polyA_pos))
    
    # add value to the full table
    fastq_table_full$polyA_final_length[match(fastq_table$name, fastq_table_full$name)] <- fastq_table$polyA_final_length
    fastq_table_full$polyA_final_pos[match(fastq_table$name, fastq_table_full$name)] <- fastq_table$polyA_final_pos
    
    # decrease polyA length for pattern matching by 1
    polyA_min_length <- polyA_min_length - 1 
    
    cat(polyA_min_length, "\n")
    
  }
  
  return(fastq_table_full)
  
}) %>% 
  dplyr::bind_rows(.)

# remove reads which have short polyA pattern or pattern not close enough (100 bp) to beggining/end of read
all_reads_polyA_filtered <- 
  all_reads_polyA %>% 
  dplyr::mutate(seq_length = nchar(fastq_seq), 
                polyA_relative_pos = seq_length - polyA_final_pos - polyA_final_length) %>% 
  dplyr::select(-fastq_seq) %>% 
  dplyr::filter((polyA_relative_pos < 100) | (polyA_final_pos < 100)) %>% 
  dplyr::filter(polyA_final_length > 10)
  
### plot as jitterplot + boxplot
ggplot(data = all_reads_polyA_filtered, aes(fastq_file, polyA_final_length)) +
  geom_boxplot(show.legend = F, outlier.shape = NA) +
  geom_jitter(aes(colour = fastq_file), size = 0.5, width = 0.25, show.legend = F) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("sample") +
  ylab("polyA length") +
  ggsave(filename = file.path(outpath, "polyA_length_all_reads_MII_GV.pdf"), width = 25, height = 15)
