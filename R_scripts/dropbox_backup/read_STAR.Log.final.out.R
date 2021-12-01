### INFO: reads and cleans STAR Log.final.out mapping statistics table
### DATE: 05. 10. 2017.
### AUTHOR: Filip Horvat

######################################################## LIBRARIES
library(dplyr)
library(stringr)
library(magrittr)
library(readr)
library(tidyr)

######################################################## FUNCTIONS
read_STAR.Log.final.out <- function(path, reshape = FALSE){
  
  # read and clean table
  STAR_df <- 
    suppressWarnings(readr::read_delim(file = path, delim = "\t", col_names = c("stat", "value"), skip = 5, trim_ws = TRUE)) %>% 
    dplyr::filter(!is.na(value)) %>%                                                  
    dplyr::mutate(stat = stringr::str_replace(stat, " \\|$", ""),                     
                  stat = stringr::str_replace(stat, "Number of reads m", "M"),        
                  stat = stringr::str_replace(stat, " reads number", ""),             
                  value = stringr::str_replace(value, "%", ""),                       
                  value = as.numeric(value))
  
  # calculate number of unmapped reads
  columns_mapped <- c("Uniquely mapped",  "Mapped to multiple loci")
  STAR_df %<>% 
    dplyr::summarize(value = value[stat == "Number of input reads"] - sum(value[stat %in% columns_mapped]),
                     stat = "Unmapped reads") %>%
    dplyr::bind_rows(STAR_df) %>%
    dplyr::select(stat, value)
  
  # reshape 
  if(reshape){
    
    columns_out <- c(columns_mapped,  "Unmapped reads", "Number of input reads")
    STAR_df %<>% 
      dplyr::filter(stat %in% columns_out) %>%
      dplyr::mutate(stat = factor(stat, levels = columns_out)) %>%
      tidyr::spread(stat, value)  %>%
      dplyr::mutate_all(funs(as.integer)) %>% 
      dplyr::mutate(log_id = path %>% basename(.) %>% stringr::str_replace(., ".Log.final.out", "")) %>% 
      magrittr::set_colnames(colnames(.) %>%  tolower(.) %>% str_replace_all(., " ", "_")) %>% 
      dplyr::select(ncol(.), dplyr::everything())   
 
  }
  
  return(STAR_df)
  
}

