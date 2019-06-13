### INFO: .Rprofile defining few usefull functions and options
### DATE: Tue Sep 18 14:03:03 2018
### AUTHOR: Filip Horvat

######################################################## WORKING DIRECTORY

######################################################## LIBRARIES

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES

######################################################## READ DATA

######################################################## MAIN CODE
### set new hidden enviorment
.env <- new.env()

### add functions to new enviorment
# sets screen to max. width
.env$wideScreen <- function(howWide = Sys.getenv("COLUMNS")) {
  options(width = as.integer(howWide))
}

# full head of tibble or data.table
.env$headt <- function(df, n = 5) {
  head(as.data.frame(df), n = n)
}

# head and tails
.env$ht <- function(df, n = 5) {
  as.data.frame(rbind(head(df, n = n), tail(df, n = n)))
}

# converts to data.frame
.env$asdf <- function(df) {
  as.data.frame(df)
}

# prints vector in one column
.env$prnt <- function(x) {
  cat(x, sep = "\n")
}

# parses arguments from command line into named vector
.env$parseCommandLineArguments <- function(args) {
  
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

### attaches hidden enviorment
attach(.env)

### fortunes on start
if(interactive())
  try(fortunes::fortune(), silent = TRUE)
