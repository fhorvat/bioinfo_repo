#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: get expression in ENCODE mouse dataset
### DATE: Sun Jun 24 16:14:35 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/repeat_expression.20180730")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(xlsx)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1/vfranke"

# get path of repeatMasker
rmsk_path <- list.files(path = genome_dir, pattern = "rmsk\\..*clean.*.gz$", full.names = T, recursive = T)

# mapped path
bam_path <- list.files(path = "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mesAur1_vfranke", pattern = ".*total.bam$", full.names = T)

######################################################## READ DATA
# read repeatMasker
rmsk_gr <-
  read_delim(file = rmsk_path, delim = "\t", col_types = cols(start = col_double(), end = col_double())) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

######################################################## MAIN CODE
### get count of reads, save summarizedExperiment as RDS
# load bam file list in memory
bamfile_chunk <- Rsamtools::BamFileList(bam_path, yieldSize = 2000000)

# register workers for parallel counting
BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))

# summarize overlaps
se <- GenomicAlignments::summarizeOverlaps(features = rmsk_gr,
                                           reads = bamfile_chunk,
                                           mode = "Union",
                                           singleEnd = FALSE,
                                           ignore.strand = TRUE)

# save summarizeOverlaps
saveRDS(se, file = file.path(outpath, "golden_hamster.GV.repeatMasker.summarizeExperiment.RDS"))

# load summarizeOverlaps
se <- readRDS(file = file.path(outpath, "golden_hamster.GV.repeatMasker.summarizeExperiment.RDS"))

### summarize count per repeat family
# get repeatMasker data.frame
rmsk_df <- 
  rmsk_gr %>% 
  as.data.frame(.) %>% 
  tibble::as.tibble(.) %>% 
  dplyr::select(-width)

# get name-family-class relations
rmsk_relations <- 
  rmsk_df %>% 
  dplyr::select(repName, repFamily, repClass) %>% 
  distinct(repName, .keep_all = T)

# get number of insertions per repName
rmsk_repName_insertions <- 
  rmsk_df %>% 
  dplyr::group_by(repName) %>% 
  dplyr::summarise(insertion_num = n())

# get summarizeExperiment as data.frame
se_df <- 
  assay(se) %>% 
  as.data.frame(.) %>% 
  as.tibble(.) %>% 
  magrittr::set_colnames(., str_remove_all(colnames(.), ".PE.total.bam")) %>% 
  dplyr::transmute(hamster_GV_count = round(((s_GV_Hamster_r1 + s_GV_Hamster_r2) / 2 ), 2))

# join with expression counts, summarize per repName
rmsk_repName_count <- 
  dplyr::bind_cols(rmsk_df, se_df) %>% 
  dplyr::group_by(repName) %>% 
  dplyr::summarize(hamster_GV_count = sum(hamster_GV_count) %>% as.numeric(.)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., rmsk_relations, by = "repName") %>% 
  dplyr::left_join(., rmsk_repName_insertions, by = "repName") %>% 
  dplyr::select(repClass, repFamily, repName, insertion_num, hamster_GV_count) %>% 
  dplyr::arrange(desc(hamster_GV_count)) 

# summarize per repFamily
rmsk_repFamily_count <- 
  rmsk_repName_count %>% 
  dplyr::mutate(repFamily = ifelse(is.na(repFamily), repClass, repFamily)) %>% 
  dplyr::group_by(repFamily) %>% 
  dplyr::summarize(hamster_GV_count = sum(hamster_GV_count) %>% as.numeric(.), 
                   insertion_num = sum(insertion_num) %>% as.numeric(.), 
                   repClass = str_c(unique(repClass), collapse = "|")) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::select(repClass, repFamily, insertion_num, hamster_GV_count) %>% 
  dplyr::arrange(desc(hamster_GV_count)) 

# write to separate sheets in xlsx 
write.xlsx(x = rmsk_repName_count %>% as.data.frame(.), 
           file = file.path(outpath, "golden_hamster.GV.repeatMasker.counts.xlsx"), 
           sheetName = "repName_sum", 
           row.names = FALSE)
write.xlsx(x = rmsk_repFamily_count %>% as.data.frame(.) , 
           file = file.path(outpath, "golden_hamster.GV.repeatMasker.counts.xlsx"), 
           sheetName = "repFamily_sum", 
           append = TRUE, 
           row.names = FALSE)
