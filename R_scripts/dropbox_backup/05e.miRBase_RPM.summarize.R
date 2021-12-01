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
    "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/ES_DcrTrans_2012/Data/Mapped/Shrimp_mm10") %>% 
  magrittr::set_names(., str_remove(., "Data.*$") %>% basename(.))

# list of experiments
experiment_list <- names(experiment_paths)

# loop through experiments
for(experiment in experiment_list){
  
  experiment <- "T3T_DcrTrans_2011"
  
  # set summarizedOverlaps .RDS path
  se_path <- file.path(outpath, str_c("miRbase.", experiment, ".se.RDS"))
  
  # read summarizeOverlaps .RDS
  se <- readRDS(file = se_path)
  
  # library size path
  library_size_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/MosIR_expression", 
                                 str_c("mosIR.", experiment, ".counts_summary.xlsx"))
  
  # read library size df
  library_size_df <- 
    xlsx::read.xlsx(file = library_size_path, sheetName = "library_sizes") %>% 
    tibble::as.tibble(.) %>% 
    tidyr::gather(library_type, library_size, -sample_id) %>% 
    tidyr::separate(library_type, into = c("mis", "library_type"), sep = "\\.") %>% 
    tidyr::unite(sample_id, sample_id, mis, sep = ".") %>% 
    tidyr::spread(library_type, library_size)
  
  # get data.frame of counts
  mirna_rpm <- 
    assay(se) %>% 
    as.tibble(.) %>% 
    magrittr::set_colnames(str_remove(colnames(.), ".bam")) %>% 
    dplyr::mutate(mirna_id = mcols(mirna_gr)$unique_name) %>% 
    dplyr::select(mirna_id, everything()) %>% 
    dplyr::filter(!duplicated(mirna_id)) %>% 
    tidyr::gather(sample_id, count, -mirna_id) %>% 
    dplyr::filter(str_detect(sample_id, "mis_0")) %>%
    dplyr::left_join(., library_size_df, by = "sample_id") %>%
    dplyr::mutate(library_size = (`21to23nt_reads` / 1e6), 
                  rpm = count / library_size) %>% 
    dplyr::select(mirna_id, sample_id, rpm) %>% 
    tidyr::spread(sample_id, rpm) %>%
    dplyr::left_join(., mirna_df %>% dplyr::select(-strand), by = c("mirna_id" = "unique_name")) %>% 
    dplyr::select(mirna_id, coordinates, everything())
  
  # # write to separate sheets in xlsx 
  # write.xlsx(x = mirna_counts %>% 
  #              dplyr::select_at(vars("mirna_id", "coordinates", starts_with("mis_0"))) %>% 
  #              magrittr::set_colnames(., str_remove(colnames(.), "mis_[0-2]{1}_")) %>% 
  #              dplyr::arrange(desc(rowSums(.[, -c(1, 2)]))) %>% 
  #              as.data.frame(.), 
  #            file = file.path(outpath, str_c("mature_miRNA.", experiment, ".counts.xlsx")), 
  #            sheetName = "allowed_0_mismatches", 
  #            row.names = FALSE)
  # 
  # write.xlsx(x = mirna_counts %>% 
  #              dplyr::select_at(vars("mirna_id", "coordinates", starts_with("mis_1"))) %>% 
  #              magrittr::set_colnames(., str_remove(colnames(.), "mis_[0-2]{1}_")) %>% 
  #              dplyr::arrange(desc(rowSums(.[, -c(1, 2)]))) %>%
  #              as.data.frame(.), 
  #            file = file.path(outpath, str_c("mature_miRNA.", experiment, ".counts.xlsx")), 
  #            sheetName = "allowed_1_mismatches", 
  #            append = TRUE, 
  #            row.names = FALSE)
  # 
  # write.xlsx(x = mirna_counts %>% 
  #              dplyr::select_at(vars("mirna_id", "coordinates", starts_with("mis_2"))) %>% 
  #              magrittr::set_colnames(., str_remove(colnames(.), "mis_[0-2]{1}_")) %>% 
  #              dplyr::arrange(desc(rowSums(.[, -c(1, 2)]))) %>% 
  #              as.data.frame(.), 
  #            file = file.path(outpath, str_c("mature_miRNA.", experiment, ".counts.xlsx")), 
  #            sheetName = "allowed_2_mismatches", 
  #            append = TRUE, 
  #            row.names = FALSE)
  
}


