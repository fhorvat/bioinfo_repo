### INFO: 
### DATE: Tue Feb 19 19:35:38 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

######################################################## READ DATA

######################################################## MAIN CODE
### product of 2 n-digit numbers 
# reverse integer using stringi::stri_reverse function
stringiRev <- function(n){
  
  # smallest and largest n-digit number
  smallest_n <- str_c(1, rep(0, n - 1) %>% str_c(., collapse = "")) %>% as.integer(.)
  largest_n <- rep(9, n) %>% str_c(., collapse = "") %>% as.integer(.)
  
  # get products
  expand.grid(smallest_n:largest_n, smallest_n:largest_n) %>%
    as.tibble(.) %>% 
    mutate(product = (Var1 * Var2), 
           rev_product = as.integer(stringi::stri_reverse(product))) %>% 
    filter(product == rev_product) %$%
    product %>% 
    as.integer(.) %>% 
    .[which.max(.)]
  
}

#### smarter and faster way
# product 
n <- 4

# smallest and largest n-digit number
smallest_n <- str_c(1, rep(0, n - 1) %>% str_c(., collapse = "")) %>% as.integer(.)
largest_n <- rep(9, n) %>% str_c(., collapse = "") %>% as.integer(.)

### loop thorugh combinations of numbers from largest to smallest, stop when you find first palindrome in product
for (i in largest_n:smallest_n) {
  
  for (j in largest_n:smallest_n) {
    
    # get integer as string
    word <- as.character(i * j)
    
    # reverse string
    reverse <- stringi::stri_reverse(word)
    
    # palindrome check
    if(word == reverse) break
    
  }
  
  # palindrome check
  if(word == reverse) break
  
}

answer <- i * j

