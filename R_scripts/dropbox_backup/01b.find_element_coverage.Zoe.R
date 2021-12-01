### INFO: 
### DATE: Sat Oct 05 13:47:57 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/element_coverage")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# documentation path
documentation_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/annotate_LINE1"

# LINE1s with ORF info path
line1_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation/L1s_nested_ours_20180516.ZJM.tidy.csv"

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/datasets/Dicer_Mili_KO/Mapped/perfect_alignments.all_multimappers/4_merged_replicates"

# bigWig path
coverage_path <- 
  list.files(mapped_path, pattern = ".*\\.bw", full.names = T) %>% 
  .[!str_detect(., "scaled")]

# original table path
original_line1_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation/L1s_nested_ours_20180517_PS_ZJM.xlsx"

######################################################## READ DATA
# read LINE1 table
line1_tb <- 
  readr::read_csv(line1_path) %>% 
  dplyr::rename(rmsk_id = id)

# read original Zoe's table
line1_original_table <- 
  openxlsx::read.xlsx(original_line1_path, sheet = 2) %>% 
  set_colnames(., make.unique(colnames(.))) %>% 
  as_tibble(.)

# read coverage from bigWig
coverage_list <- purrr::map(coverage_path, function(path){
  
  # read bigWig file
  coverage <- 
    rtracklayer::import(path) %>% 
    .[mcols(.)$score > 0] %>% 
    reduce(., ignore.strand = T)
  
}) %>% 
  set_names(., coverage_path %>% basename(.) %>% str_remove_all(., "^s_GV_|\\.bw$"))

######################################################## MAIN CODE
### prepare 
# create GRanges
line1_gr <- 
  line1_tb %>% 
  dplyr::mutate(strand = ifelse(!(strand %in% c("+", "-")), "*", strand)) %>% 
  GRanges(.)

# get width in table
line1_tb <- as_tibble(line1_gr)

### find coverage for each LINE1 element
line1_coverage_list <- purrr::map(names(coverage_list), function(sample_id){
  
  # get individual sample coverage
  coverage <- coverage_list[[sample_id]]
  
  # find overlaps between two GRanges
  hits <- findOverlaps(line1_gr, coverage, ignore.strand = T)
  
  # extract all overalaping features from subject as list
  line1_coverage <- extractList(coverage, as(hits, "List"))
  
  # set names to diff. list, unlist
  names(line1_coverage) <- mcols(line1_gr)$rmsk_id 
  
  
  ### find longest stretch and percentage of coverage for each LINE1
  # unlist
  line1_coverage_tb <- 
    line1_coverage %>% 
    unlist(.)
  mcols(line1_coverage_tb)$rmsk_id <- names(line1_coverage_tb)
  
  # convert to data.table, summarize 
  line1_coverage_sum <- 
    line1_coverage_tb %>% 
    as.data.table(.) %>% 
    .[,.(coverage_total = sum(width), coverage_longest = max(width)), by = rmsk_id] %>% 
    as_tibble(.) %>% 
    dplyr::left_join(., line1_tb %>% dplyr::select(rmsk_id, width), by = "rmsk_id") %>% 
    dplyr::mutate(coverage_ratio = round((coverage_total / width), 3)) %>% 
    dplyr::select(rmsk_id, coverage_ratio, longest_stretch = coverage_longest) %>% 
    set_colnames(., c("rmsk_id", str_c(c("coverage_ratio", "longest_stretch"), sample_id, sep = ".")))
  
  # return 
  return(line1_coverage_sum)
  
})

# join to one table
line1_coverage_tb <- 
  purrr::reduce(c(list(line1_tb), line1_coverage_list), left_join, by = "rmsk_id") %>% 
  dplyr::mutate_at(vars(starts_with("coverage_ratio")), list(~ifelse(is.na(.), 0, .))) %>% 
  dplyr::mutate_at(vars(starts_with("longest_")), list(~ifelse(is.na(.), 0, .)))

# filter table and save
line1_out <- 
  line1_coverage_tb %>% 
  # dplyr::filter(width >= 5000, width <= 7000, 
  #               insertion_class %in% c("whole", "within")) %>% 
  dplyr::rename(id = rmsk_id) %>% 
  dplyr::select(-c(repClass, repFamily)) %T>%
  readr::write_csv(., path = file.path(outpath, "LINE1.Zoe.Dicer_Mili_KO.coverage.csv"))


### original table
# join with original table and save
line1_original_table_rpm <- 
  left_join(line1_original_table, 
            line1_out %>% dplyr::select_at(vars("id", starts_with("coverage_ratio"), starts_with("longest_stretch"))), 
            by = c("uniqName" = "id")) %>% 
  write_csv(., file.path(outpath, str_c("L1s_nested_ours_20180517_PS_ZJM.original.", "Dicer_Mili_KO", ".coverage.csv")))
