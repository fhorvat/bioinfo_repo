### INFO: 
### DATE: Thu Jul 30 18:45:22 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LTRs/expression/hamster_oocyte_Mov10l.RNAseq")

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

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# mapped path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq"
documentation_path <- file.path(base_path, "Data/Documentation")
sample_table_path <- list.files(documentation_path, "\\.sampleTable\\.csv", full.names = T)
mapped_path <- file.path(base_path, "Data/Mapped/STAR_Siomi.multimappers")
library_size_path <- file.path(mapped_path, "3_logs", "log.read_stats.txt")

# chosen LTR classes path
classes_tb_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LTRs/substitution_rate"
classes_tb_path <- file.path(classes_tb_path, "all_LTR_classes 200730.xlsx")

# expression files path
expression_path <- file.path(inpath, "expression_sum.RDS_files")
expression_path <- list.files(expression_path, ".*\\.expression\\.RDS$", full.names = T)

######################################################## READ DATA
# read sample table
sample_tb <- readr::read_csv(sample_table_path) 

# read library sizes
library_size_tb <- 
  readr::read_delim(library_size_path, delim = "\t") %>% 
  dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA)

# read table with chosen LTR classes
classes_tb <- openxlsx::read.xlsx(classes_tb_path)

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
### prepare files
# select relevant columns in classes table
classes_tb %<>% 
  as_tibble(.) %>% 
  dplyr::select(repName, repFamily, category_I, type) %>% 
  dplyr::mutate(repFamily = replace(repFamily, is.na(repFamily), "other"))

### get expression 
# add category to expression table, join with library size, get RPM, join with sample table, summarize per genotype
rpm_tb <- 
  expression_tb %>% 
  dplyr::left_join(., classes_tb, by = "repName") %>% 
  dplyr::group_by(sample_id, category_I) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., library_size_tb, by = "sample_id") %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt$"), 
                library_size = library_size / 1e6,
                RPM = count / library_size) %>%
  dplyr::select(sample_id, category_I, RPM)

# filter table
rpm_filt <- 
  rpm_tb %>% 
  dplyr::left_join(., sample_tb, by = "sample_id") %>% 
  dplyr::group_by(genotype, category_I) %>% 
  dplyr::summarise(RPM = round(mean(RPM), 3)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(!str_detect(genotype, "HET")) %>% 
  tidyr::pivot_wider(id_cols = c("genotype", "category_I"), names_from = genotype, values_from = RPM, names_prefix = "RPM.") %>% 
  dplyr::mutate_all(.funs = ~(replace(., is.na(.), 0)))


### save 
# open workbook and sheets
wb <- createWorkbook("RPM")

# add sheet
addWorksheet(wb, sheetName = "RPM")

# write to workbook
writeDataTable(wb = wb, 
               sheet = "RPM", 
               x = rpm_filt)

# save workbook to disk
saveWorkbook(wb, file.path(outpath, str_c("all_LTR_classes.200730", "RNA.oocyte", "mean_RPM", "xlsx", sep = ".")), overwrite = TRUE)
