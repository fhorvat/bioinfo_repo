### INFO: scraps http://os-akovacica-mgorica.skole.hr/ front page for titles and subtitles along with dates
### DATE: Thu Aug 22 16:18:52 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/other/web_scraping")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(purrr)

library(rvest)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

######################################################## READ DATA

######################################################## MAIN CODE
# set url
url <- "http://os-akovacica-mgorica.skole.hr"

# read .html from the website
webpage <- read_html(url)

# set all .CSS selectors you want to scrape
selectors <- c(".naslov", ".autor")

### get all selectors and join to one table
page_tb <-
  purrr::map(selectors, function(selector){
    
    # get all titles
    data_html <- html_nodes(webpage, selector)
    
    # convert to text
    data_txt <- html_text(data_html)
    
    # return
    return(data_txt)
    
  }) %>% 
  dplyr::bind_cols(.) %>% 
  set_colnames(., c("naslov", "datum"))

# tidy table
page_tidy <- 
  page_tb %>% 
  dplyr::mutate(datum =  str_remove_all(datum, "^.* datum: |\t|\n") %>% as.Date(., format = "%d. %m. %Y. %H:%M")) %>% 
  dplyr::filter(datum >= "2018-09-09") %>% 
  arrange(datum) %>% 
  dplyr::mutate(datum = format(datum, "%d. %m. %Y.")) %>% 
  dplyr::mutate(output = str_c(datum, " - ", naslov)) %$%
  output %T>% 
  readr::write_lines(., file.path(outpath, "naslovi_2019_2020.20200505.txt"))

