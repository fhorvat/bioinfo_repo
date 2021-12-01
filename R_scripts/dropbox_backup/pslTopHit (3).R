#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: takes psl and scored psl as input, outputs top hit psl
### DATE: Tue Jul 23 12:26:22 2019
### AUTHOR: Filip Horvat

######################################################## WORKING DIRECTORY

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# blat .psl path
blat_psl_path <- args$psl

# blat .psl path with score
blat_psl_score_path <- str_replace(blat_psl_path, "\\.psl$", ".score.psl")

######################################################## READ DATA
# read .psl
psl <- readr::read_delim(blat_psl_path, delim = "\t", skip = 5, col_types = cols(.default = "c"),
                         col_names = c("match", "mismatch", "rep_match", "Ns", 
                                       "Q_gap_count", "Q_gap_bases", "T_gap_count", "T_gap_bases", "strand", 
                                       "qName", "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd",
                                       "blockCount", "blockSizes", "qStarts", "tStarts"))

# read .psl header
psl_header <- readr::read_lines(blat_psl_path, n_max = 5)

# read .psl with score
psl_score <- readr::read_delim(file = blat_psl_score_path, delim = "\t", col_types = cols(.default = "c"),
                               col_names = c("tName", "tStart", "tEnd", "q_coordinates", "score", "percentIdentity"))

######################################################## MAIN CODE
# if the .psl file isn't empty get top score and write
if(nrow(psl) > 0) {
  
  # join score with psl, get top 1 hit by score, collapse
  psl_top <- 
    psl %>% 
    left_join(., psl_score %>% dplyr::select(-c(q_coordinates, percentIdentity)), by = c("tName", "tStart", "tEnd")) %>% 
    dplyr::mutate(score = as.numeric(score)) %>% 
    top_n(1, score) %>% 
    dplyr::select(-score) %>% 
    dplyr::slice(1) %>% 
    unlist(., use.names = FALSE) %>% 
    str_c(., collapse = "\t")
  
  # write along with header
  readr::write_lines(c(psl_header, psl_top), path = str_replace(blat_psl_path, "\\.all\\.psl$", ".psl"))
  
}else{
  
  # write just header
  readr::write_lines(c(psl_header), path = str_replace(blat_psl_path, "\\.all\\.psl$", ".psl"))
  
}
