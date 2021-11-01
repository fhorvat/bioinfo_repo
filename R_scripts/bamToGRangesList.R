### INFO: reads .bam file and outputs GenomicRangesList (ranges are with respect to CIGAR string)
### DATE: 04. 10. 2017.
### AUTHOR: Filip Horvat

######################################################## FUNCTIONS
library(magrittr)
library(GenomicAlignments)
library(GenomicRanges)

######################################################## FUNCTIONS
bamToGRangesList <- function(bam_path){
  
  GenomicAlignments::readGAlignmentsList(bam_path, use.names = TRUE) %>% 
    unlist(.) %>% 
    GenomicRanges::grglist(.)
  
}
