### INFO: .Rprofile defining few usefull functions and options
### DATE: Tue Sep 18 14:03:03 2018
### AUTHOR: Filip Horvat

######################################################## WORKING DIRECTORY

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES

######################################################## READ DATA

######################################################## MAIN CODE
### set new hidden enviorment
.env <- new.env()

### adds variable called "pwd" to new enviorment
.env$pwd <- getwd()

### add functions to new enviorment
## sets screen to max. width
.env$wideScreen <- function(howWide = Sys.getenv("COLUMNS")) {
  options(width = as.integer(howWide))
}

## full head of tibble or data.table
.env$headt <- function(df, n = 5) {
  head(as.data.frame(df), n = n)
}

## head and tails
.env$ht <- function(df, n = 5) {
  as.data.frame(rbind(head(df, n = n), tail(df, n = n)))
}

## converts to data.frame
.env$asdf <- function(df) {
  as.data.frame(df)
}

## prints vector in one column
.env$prnt <- function(x) {
  cat(x, sep = "\n")
}


### INFO: parses arguments from command line into named vector
### DATE: 17. 06. 2019.
### AUTHOR: Filip Horvat
.env$parseCommandLineArguments <- function(args) {
  
  # libraries
  require(stringr)
  require(magrittr)
  require(purrr)

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


### INFO: expands GRanges up- and down-stream
### DATE: 2018-05-30
### AUTHOR: Devon Ryan, https://bioinformatics.stackexchange.com/questions/4390/expand-granges-object-different-amounts-upstream-vs-downstream
.env$expandRange <- function(x, upstream=2000, downstream=1000) {
  
  strand_is_minus <- strand(x) == "-"
  on_plus <-- which(!strand_is_minus)
  on_minus <- which(strand_is_minus)
  
  start(x)[on_plus] <- start(x)[on_plus] - upstream
  start(x)[on_minus] <- start(x)[on_minus] - downstream
  
  end(x)[on_plus] <- end(x)[on_plus] + downstream
  end(x)[on_minus] <- end(x)[on_minus] + upstream
  
 return(x)

}

### attaches hidden enviorment
attach(.env)
