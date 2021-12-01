### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/siRNA_editing")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# MosIR reads path
mosir_reads_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Mapped/STAR_mm10.pCAG_EGFP_MosIR.new/7_MosIR_reads"

# MosIR tsv
mosir_tsv_path_list <- 
  list.files(mosir_reads_path, pattern = ".*\\.tsv$", full.names = T) %>% 
  set_names(., str_remove(basename(.), ".SE.18to30nt.MosIR.tsv$"))

# MosIR bam
mosir_bam_path_list <- 
  list.files(mosir_reads_path, pattern = ".*\\.bam$", full.names = T) %>% 
  set_names(., str_remove(basename(.), ".SE.18to30nt.MosIR.bam$"))
  
######################################################## READ DATA

######################################################## MAIN CODE
# make widths dummy table
dummy_widths <- 
  expand.grid(width = 18:30, 
              base = c("A", "T", "C", "G"), 
              ref = c("A", "T", "C", "G"), 
              stringsAsFactors = F) %>% 
  as_tibble(.) 

### get table of mismatches of reads mapped to MosIR
siRNA_editing <- purrr::map(names(mosir_tsv_path_list), function(sample_name){
  
  # get paths to files
  mosir_tsv_path <- mosir_tsv_path_list[sample_name]
  mosir_bam_path <- mosir_bam_path_list[sample_name]
  
  # get read mismatches
  mosir_read_mismatches <- 
    suppressMessages(suppressWarnings(read_delim(file = mosir_tsv_path, delim = "\t"))) %>% 
    dplyr::rename(READ_NAME = `#READ_NAME`) %>% 
    magrittr::set_colnames(., tolower(colnames(.))) %>% 
    dplyr::filter(!str_detect(read_name, "^:.*"), 
                  base != ref, 
                  op != "S")
  
  # get read widths
  mosir_read_widths <- 
    GenomicAlignments::readGAlignmentsList(file = mosir_bam_path, 
                                           use.names = TRUE, 
                                           param = ScanBamParam(flag = scanBamFlag(isMinusStrand = F))) %>% 
    unlist(.) %>% 
    as.data.frame(.) %>% 
    as_tibble(., rownames = "read_name") %>% 
    dplyr::select(read_name, width)
    
  # join tables, spread
  mosir_editing <- 
    mosir_read_mismatches %>% 
    dplyr::inner_join(., mosir_read_widths, by = "read_name") %>% 
    dplyr::group_by(width) %>% 
    dplyr::count(base, ref) %>% 
    dplyr::right_join(., dummy_widths) %>% 
    dplyr::arrange(width, base, ref) %>% 
    dplyr::mutate(n = replace(n, is.na(n), 0)) %>% 
    tidyr::spread(ref, n) %>% 
    dplyr::mutate(sample_name = sample_name) %>% 
    dplyr::select(sample_name, everything(.))
  
})

# join to one table, write
dplyr::bind_cols(siRNA_editing) %>%
  openxlsx::write.xlsx(., file = file.path(outpath, str_c("Eliska_mESC_MosIR", "MosIR_siRNA_editing.xlsx")))

