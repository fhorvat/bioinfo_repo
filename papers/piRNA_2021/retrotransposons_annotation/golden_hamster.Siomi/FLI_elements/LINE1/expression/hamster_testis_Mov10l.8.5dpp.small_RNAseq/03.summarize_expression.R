### INFO: 
### DATE: Thu Jul 30 18:45:22 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/expression/hamster_testis_Mov10l.8.5dpp.small_RNAseq")

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

library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# mapped path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.small_RNAseq"
documentation_path <- file.path(base_path, "Data/Documentation")
sample_table_path <- list.files(documentation_path, "\\.sampleTable\\.csv", full.names = T)
mapped_path <- file.path(base_path, "Data/Mapped/STAR_Siomi.multimappers")
library_size_path <- file.path(mapped_path, "library_sizes.txt")

# expression files path
expression_path <- file.path(inpath, "expression_sum.RDS_files")
expression_path <- list.files(expression_path, ".*\\.expression\\.RDS$", full.names = T)

######################################################## READ DATA
# read sample table
sample_tb <- readr::read_csv(sample_table_path) 

# read library sizes
library_size_tb <- readr::read_delim(library_size_path, delim = "\t", col_names = c("sample_id", "library_size"))

# read expression values
expression_tb <- purrr::map(expression_path, function(path){
  
  # get sample name
  sample_name <- path %>% basename(.) %>% str_remove(., "\\.expression\\.RDS")
  
  # read table
  tmp_tb <- 
    path %>% 
    readRDS(.) %>% 
    dplyr::mutate(sample_id = sample_name)
  
  # return
  return(tmp_tb)
  
}) %>% 
  dplyr::bind_rows(.)

######################################################## MAIN CODE
### set element name
element_name <- "IAP"

### get expression 
# join with library size, get RPM
rpm_tb <- 
  expression_tb %>% 
  dplyr::left_join(., library_size_tb, by = "sample_id") %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt$"), 
                library_size = library_size / 1e6,
                RPM = count / library_size) %>%
  dplyr::select(sample_id, sense, read_width, RPM)

### sum different read lengths
# read lengths list
read_lengths <- list("19to32nt" = c(19, 32), "19to24nt" = c(19, 24), "25to32nt" = c(25, 32), "24to31nt" = c(24, 31))

# open workbook and sheets
wb <- createWorkbook("RPM")

# loop through read lengths list, save as different sheets in .xlsx
purrr::map(names(read_lengths), function(rlen){
  
  # join with sample table, summarize per genotype/read length
  rpm_mean <- 
    rpm_tb %>% 
    dplyr::filter(read_width >= min(read_lengths[[rlen]]), read_width <= max(read_lengths[[rlen]])) %>% 
    dplyr::group_by(sample_id, sense) %>%
    dplyr::summarise(RPM = sum(RPM)) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::left_join(., sample_tb, by = "sample_id") %>% 
    tidyr::unite("genotype_age", genotype, age, sep = "_") %>% 
    dplyr::group_by(genotype_age, sense) %>% 
    dplyr::summarise(RPM = round(mean(RPM), 3)) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::filter(!str_detect(genotype_age, "HET")) %>% 
    tidyr::pivot_wider(id_cols = "genotype_age", names_from = sense, names_prefix = "RPM.", values_from = RPM) %>% 
    dplyr::mutate_all(.funs = ~(replace(., is.na(.), 0)))
  
  # add sheet
  addWorksheet(wb, sheetName = rlen)
  
  # write to workbook
  writeDataTable(wb = wb, 
                 sheet = rlen, 
                 x = rpm_mean)
  
  # save workbook to disk
  saveWorkbook(wb, file.path(outpath, str_c(element_name, "FLI_elements.RPM", "smallRNA", "testis", "8.5dpp", "xlsx", sep = ".")), overwrite = TRUE)
  
  # return 
  return(rlen)
  
})


