### INFO: reads and parses BLAST format 6 output
### DATE: 14. 09. 2017.
### AUTHOR: Filip Horvat

parseBLAST6tbl <- function(blast_tbl_path){
  
  # read and parse table
  blast6 <- 
    readr::read_delim(file = blast_tbl_path, delim = "\t", 
                      col_names = c("query_id", "subject_id", "identity_perc", "alignment_length", 
                                    "mismatches", "gap", "openings", 
                                    "query_start", "query_end", "subject_start", "subject_end", 
                                    "e_value", "bit_score")) %>%
    dplyr::mutate(query_fullName = readr::str_c(seqnames, ":", start, "-", end, ":", strand))
  
}
