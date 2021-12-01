### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/class_reads")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(xlsx)

library(GenomicAlignments)
library(GenomicRanges)
library(GenomicFeatures)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# path to plasmid coordinates
plasmid_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/mm10.all_plasmids/plasmid_sequences/all_plasmids.bed"

######################################################## READ DATA
# read plasmid bed
plasmid_df <- readr::read_delim(plasmid_path, delim = "\t", col_names = c("seqnames", "start", "end"))

######################################################## MAIN CODE
# plasmid names
plasmid_names <- plasmid_df$seqnames

# set categories order
read_classes <- c("IR", "miRNA.mature.sense", "miRNA.other.sense", "protein_coding.sense", 
                  "rRNA", "SINE", "LINE", "LTR", "other_repeat", "annotated_pseudogene", 
                  "other", "not_annotated")

### adds different library sizes (whole mapped library, 18 to 30 nt reads, 21 to 23 nt reads) to results tables
# experiment paths
experiment_paths <- 
  "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/NIH3T3_transfected.2018/Data/Mapped/STAR_mm10.plasmids" %>% 
  magrittr::set_names(., str_remove(., "Data.*$") %>% basename(.))

# list of experiments
experiment_list <- names(experiment_paths)

# loop through experiments
experiment <- experiment_list[1]

# get mapped path 
mapped_path <- experiment_paths[experiment]

# classified reads path
read_class_path <- list.files(path = file.path(mapped_path, "5_class_reads"), pattern = "read_class.*mis_[0-2]{1}.txt", full.names = T)

# read class reads table, gets different mismatch counts
read_class_df <-
  purrr::map(read_class_path, readr::read_delim, delim = "\t") %>%
  dplyr::bind_rows(.) %>% 
  dplyr::rename(alignment_count = count) %>% 
  dplyr::select(-experiment) %>% 
  tidyr::separate(sample_id, c("sample_id", "mismatch"), sep = ".mis_") %>% 
  dplyr::mutate(mismatch = str_c("mis_", mismatch), 
                read_group = replace(read_group, read_group %in% plasmid_names, "IR")) %>% 
  tidyr::spread(mismatch, alignment_count) %>% 
  dplyr::mutate_at(vars(starts_with("mis")), funs(replace(., is.na(.), 0))) %>% 
  dplyr::mutate(mis_1 = mis_1 + mis_0, 
                mis_2 = mis_2 + mis_1) %>% 
  tidyr::gather(mismatch_allowed, alignment_count, mis_0:mis_2) %>% 
  tidyr::unite(tmp, mismatch_allowed, read_group, sep = ".") %>% 
  tidyr::spread(tmp, alignment_count)

# write to separate sheets in xlsx 
write.xlsx(x = read_class_df %>%
             dplyr::select(sample_id, str_c("mis_0.", read_classes)) %>% 
             magrittr::set_colnames(., str_remove(colnames(.), "mis_[0-2]{1}.")) %>% 
             as.data.frame(.), 
           file = file.path(outpath, str_c("class_reads.hierarchy", experiment, ".hierarchy_counts.xlsx")), 
           sheetName = "allowed_0_mismatches", 
           row.names = FALSE)

write.xlsx(x = read_class_df %>%
             dplyr::select(sample_id, str_c("mis_1.", read_classes)) %>% 
             magrittr::set_colnames(., str_remove(colnames(.), "mis_[0-2]{1}.")) %>% 
             as.data.frame(.), 
           file = file.path(outpath, str_c("class_reads.hierarchy", experiment, ".hierarchy_counts.xlsx")), 
           sheetName = "allowed_1_mismatches", 
           append = TRUE, 
           row.names = FALSE)

write.xlsx(x = read_class_df %>%
             dplyr::select(sample_id, str_c("mis_2.", read_classes)) %>% 
             magrittr::set_colnames(., str_remove(colnames(.), "mis_[0-2]{1}.")) %>% 
             as.data.frame(.), 
           file = file.path(outpath, str_c("class_reads.hierarchy", experiment, "counts.xlsx")), 
           sheetName = "allowed_2_mismatches", 
           append = TRUE, 
           row.names = FALSE)


