### INFO: reads and cleans SJ.out.tab output from STAR aligner
### DATE: 14. 09. 2017.
### AUTHOR: Filip Horvat

######################################################## LIBRARIES
library(magrittr)
library(readr)
library(dplyr)
library(stringr)

######################################################## FUNCTIONS
read_SJout <- function(path){
  
  # set info about strand/intron motifs coding
  all_strands <- c("*", "+", "-") %>% magrittr::set_names(as.character(0:2))
  all_intron_motifs <- c("non-canonical", "GT/AG", "CT/AC", "GC/AG", "CT/GC", "AT/AC", "GT/AT") %>% magrittr::set_names(as.character(0:6))
  
  # read SJ.out table from STAR
  sjout <- 
    readr::read_delim(file = path, delim = "\t", col_names = c("seqnames", "start", "end", "strand", "intron_motif", "annot", "n_uniq", "n_multi", "overhang")) %>%
    dplyr::mutate(strand = all_strands[as.character(strand)], 
                  intron_motif = all_intron_motifs[as.character(intron_motif)], 
                  SJ_fullName = stringr::str_c(seqnames, ":", start, "-", end, ":", strand))

  return(sjout)

}


