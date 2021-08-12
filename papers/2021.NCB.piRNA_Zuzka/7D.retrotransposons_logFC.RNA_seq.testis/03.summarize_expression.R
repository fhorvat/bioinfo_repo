### INFO: 
### DATE: Thu Jul 30 18:45:22 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/raw_rmsk.expression/expression/small_RNAseq")

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
# set outpath
outpath <- getwd()

# chosen LTR classes path
classes_tb_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LTRs/substitution_rate"
classes_tb_path <- file.path(classes_tb_path, "all_LTR_classes 200730.xlsx")

######################################################## READ DATA
# read table with chosen LTR classes
classes_tb <- openxlsx::read.xlsx(classes_tb_path)

#################################### MAIN CODE
### prepare files
# select relevant columns in classes table
classes_tb %<>% 
  as_tibble(.) %>% 
  dplyr::select(repName, repFamily, category_I, type) %>% 
  dplyr::mutate(repFamily = replace(repFamily, is.na(repFamily), "other"))

# set list of dataset names
dataset_name_list <- c("hamster_testis_Mov10l.8.5dpp.small_RNAseq", 
                       "hamster_testis_Mov10l.small_RNAseq", 
                       "hamster_testis_Mov10l.small_RNAseq.reseq", 
                       "hamster_testis_Mov10l.8.5dpp.run_2.small_RNAseq")

# get expression values for different datasets
rpm_tb_all <- purrr::map(dataset_name_list, function(dataset_name){
  
  ### set paths
  # mapped path
  base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets", dataset_name)
  documentation_path <- file.path(base_path, "Data/Documentation")
  sample_table_path <- list.files(documentation_path, ".*\\.sampleTable.csv", full.names = T)
  mapped_path <- file.path(base_path, "Data/Mapped/STAR_Siomi.multimappers")
  library_size_path <- file.path(mapped_path, "4_library_size", "library_sizes.txt")
  
  # expression files path
  inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly"
  inpath <- file.path(inpath, "raw_rmsk.expression/expression/small_RNAseq", dataset_name, "expression_sum.RDS_files")
  expression_path <- list.files(path = inpath, pattern = ".*\\.expression\\.RDS$", full.names = T)
  
  
  ## read data
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
  
  
  ### calculate RPMs
  # add category to expression table, join with library size, get RPM, join with sample table, summarize per genotype
  rpm_tb <- 
    expression_tb %>% 
    dplyr::left_join(., classes_tb, by = "repName") %>% 
    dplyr::filter(!is.na(category_I)) %>% 
    dplyr::group_by(sample_id, category_I, sense, read_width) %>% 
    dplyr::summarise(count = sum(count)) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::left_join(., library_size_tb, by = "sample_id") %>% 
    dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt$"), 
                  library_size = library_size / 1e6,
                  RPM = count / library_size) %>%
    dplyr::select(sample_id, category_I, sense, read_width, RPM)
  
  # return
  return(rpm_tb)
  
}) %>% 
  dplyr::bind_rows(.)


### sum different read lengths
# read lengths list
read_lengths <- list("19to32nt" = c(19, 32), "19to24nt" = c(19, 24), "25to32nt" = c(25, 32), "24to31nt" = c(24, 31))

# open workbook and sheets
wb_individual <- createWorkbook("RPM")
wb_mean <- createWorkbook("RPM")

# loop through read lengths list, save as different sheets in .xlsx
purrr::map(names(read_lengths), function(rlen){
  
  # filter table
  rpm_long <- 
    rpm_tb_all %>% 
    dplyr::filter(read_width >= min(read_lengths[[rlen]]), read_width <= max(read_lengths[[rlen]])) %>% 
    dplyr::right_join(., sample_tb, by = "sample_id") %>% 
    dplyr::group_by(sample_id, category_I, sense) %>%
    dplyr::summarise(RPM = sum(RPM), 
                     age = unique(age), 
                     genotype = unique(genotype)) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::arrange(sense, age, genotype)
  
  # get wide table with individual samples RPM
  rpm_individual <- 
    rpm_long %>% 
    tidyr::pivot_wider(id_cols = category_I, names_from = c(sample_id, sense), values_from = RPM) %>% 
    dplyr::mutate_all(.funs = ~(replace(., is.na(.), 0)))
  
  # write to workbook
  addWorksheet(wb_individual, sheetName = rlen)
  writeDataTable(wb = wb_individual, 
                 sheet = rlen, 
                 x = rpm_individual)
  
  # get wide table with genotype/age mean RPM
  rpm_mean <- 
    rpm_long %>% 
    dplyr::group_by(sense, age, genotype, category_I) %>% 
    dplyr::summarise(RPM = round(mean(RPM), 3)) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::filter(str_detect(genotype, "WT")) %>% 
    tidyr::pivot_wider(id_cols = c(category_I), names_from = c(genotype, age, sense), values_from = RPM, names_sep = ".") %>% 
    dplyr::mutate_all(.funs = ~(replace(., is.na(.), 0)))
  
  # write to workbook
  addWorksheet(wb_mean, sheetName = rlen)
  writeDataTable(wb = wb_mean, 
                 sheet = rlen, 
                 x = rpm_mean)
  
  # return 
  return(rlen)
  
})

# save workbooks to disk
saveWorkbook(wb_individual, file.path(outpath, str_c("all_LTR_classes", "small_RNAseq", "all_samples_RPM", "20210104", "xlsx", sep = ".")), overwrite = TRUE)
saveWorkbook(wb_mean, file.path(outpath, str_c("all_LTR_classes", "small_RNAseq", "mean_RPM", "20210104", "xlsx", sep = ".")), overwrite = TRUE)
