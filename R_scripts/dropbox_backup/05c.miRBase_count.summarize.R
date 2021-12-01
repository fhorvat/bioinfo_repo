### INFO: counts reads over mature miRNAs (from miRbase .gff)
### DATE: Thu Aug 16 17:26:28 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/miRNA_expression")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(xlsx)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# miRbase path
mirbase_path <- list.files(genome_path, pattern = "miRBase.*gff3", full.names = T)

######################################################## READ DATA
# read miRbase gtf
mirbase_gr <- rtracklayer::import.gff(con = mirbase_path) 

######################################################## MAIN CODE
# get ranges of mature miRNA
mirna_gr <- mirbase_gr[mcols(mirbase_gr)$type == "miRNA"]
mcols(mirna_gr) <- mcols(mirna_gr)[, c("Name", "Derives_from")]
mcols(mirna_gr)$unique_name <- str_c(mcols(mirna_gr)$Name, ".", mcols(mirna_gr)$Derives_from)
  
# create data.frame with miRNA coordinates
mirna_df <- 
  mirna_gr %>% 
  as.data.frame(.) %>% 
  as.tibble(.) %>% 
  tidyr::unite(coordinates, seqnames, start, end, sep = " ") %>% 
  dplyr::select(unique_name, coordinates, strand) %>% 
  dplyr::filter(!duplicated(unique_name))
  
### read summarizedOverlaps .RDS for all experiments and writes results table
# experiment paths
experiment_paths <- 
  c("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/T3T_DcrTrans_2011/Data/Mapped/Shrimp_mm10", 
    "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/ES_DcrTrans_2012/Data/Mapped/Shrimp_mm10", 
    "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/T3T_mESC_MosIR.2016/Data/Mapped/STAR_mm10",
    "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Mapped/STAR_mm10.pCAG_EGFP_MosIR", 
    "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/NIH3T3_transfected.2018/Data/Mapped/STAR_mm10.all_plasmid") %>% 
  magrittr::set_names(., str_remove(., "Data.*$") %>% basename(.))

# list of experiments
experiment_list <- names(experiment_paths)

# loop through experiments
for(experiment in experiment_list){

  experiment <- "NIH3T3_transfected.2018"
  
  # set summarizedOverlaps .RDS path
  se_path <- file.path(outpath, str_c("miRbase.", experiment, ".se.RDS"))
  
  # read summarizeOverlaps .RDS
  se <- readRDS(file = se_path)
  
  # get data.frame of counts
  mirna_counts <- 
    assay(se) %>% 
    as.tibble(.) %>% 
    magrittr::set_colnames(str_remove(colnames(.), ".bam")) %>% 
    dplyr::mutate(mirna_id = mcols(mirna_gr)$unique_name) %>% 
    dplyr::select(mirna_id, everything()) %>% 
    dplyr::filter(!duplicated(mirna_id)) %>% 
    tidyr::gather(sample_id, count, -mirna_id) %>% 
    tidyr::separate(sample_id, c("sample_id", "mismatch"), sep = ".mis_") %>% 
    dplyr::mutate(mismatch = str_c("mis_", mismatch)) %>% 
    tidyr::spread(mismatch, count) %>% 
    dplyr::mutate_at(vars(starts_with("mis")), funs(replace(., is.na(.), 0))) %>% 
    dplyr::mutate(mis_1 = mis_1 + mis_0, 
                  mis_2 = mis_2 + mis_1) %>% 
    dplyr::arrange(sample_id) %>% 
    tidyr::gather(mismatch_allowed, count, mis_0:mis_2) %>%
    tidyr::unite(temp, mismatch_allowed, sample_id) %>%
    tidyr::spread(temp, count) %>%
    dplyr::mutate_at(vars(starts_with("mis_")), funs(replace(., is.na(.), 0))) %>% 
    dplyr::left_join(., mirna_df %>% dplyr::select(-strand), by = c("mirna_id" = "unique_name")) %>% 
    dplyr::select(mirna_id, coordinates, everything())
  
  # write to separate sheets in xlsx 
  write.xlsx(x = mirna_counts %>% 
               dplyr::select_at(vars("mirna_id", "coordinates", starts_with("mis_0"))) %>% 
               magrittr::set_colnames(., str_remove(colnames(.), "mis_[0-2]{1}_")) %>% 
               dplyr::arrange(desc(rowSums(.[, -c(1, 2)]))) %>% 
               as.data.frame(.), 
             file = file.path(outpath, str_c("mature_miRNA.", experiment, ".counts.xlsx")), 
             sheetName = "allowed_0_mismatches", 
             row.names = FALSE)
  
  write.xlsx(x = mirna_counts %>% 
               dplyr::select_at(vars("mirna_id", "coordinates", starts_with("mis_1"))) %>% 
               magrittr::set_colnames(., str_remove(colnames(.), "mis_[0-2]{1}_")) %>% 
               dplyr::arrange(desc(rowSums(.[, -c(1, 2)]))) %>%
               as.data.frame(.), 
             file = file.path(outpath, str_c("mature_miRNA.", experiment, ".counts.xlsx")), 
             sheetName = "allowed_1_mismatches", 
             append = TRUE, 
             row.names = FALSE)
  
  write.xlsx(x = mirna_counts %>% 
               dplyr::select_at(vars("mirna_id", "coordinates", starts_with("mis_2"))) %>% 
               magrittr::set_colnames(., str_remove(colnames(.), "mis_[0-2]{1}_")) %>% 
               dplyr::arrange(desc(rowSums(.[, -c(1, 2)]))) %>% 
               as.data.frame(.), 
             file = file.path(outpath, str_c("mature_miRNA.", experiment, ".counts.xlsx")), 
             sheetName = "allowed_2_mismatches", 
             append = TRUE, 
             row.names = FALSE)
  
}


