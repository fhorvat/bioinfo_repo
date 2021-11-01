### INFO: transforms arguments from command line into named vector 
### DATE: Mon May 20 16:55:44 2019
### AUTHOR: Filip Horvat

######################################################## LIBRARIES
library(magrittr)
library(stringr)
library(purrr)

######################################################## FUNCTIONS
parseCommandLineArguments <- function(args){
  
  # get named list of arguments
  args_list <-
    args %>%
    unlist(.) %>%
    stringr::str_c(., collapse = " ") %>%
    stringr::str_split(., "--") %>%
    unlist(.) %>%
    .[-1] %>% 
    purrr::map(., function(x){
      
      # split input
      input_split <- 
        x %>% 
        stringr::str_squish(.) %>% 
        stringr::str_split(., " ") %>% 
        unlist(.)
      
      # get arguments
      input_args <- 
        input_split[-1] %>% 
        list(.) %>% 
        set_names(., input_split[1])
      
      # return
      return(input_args)
      
    }) %>% 
    unlist(., recursive = F)
  
}
