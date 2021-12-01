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
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# get path of repeatMasker
rmsk_path <- list.files(path = genome_dir, pattern = "rmsk\\..*clean.*.gz$", full.names = T, recursive = T)

# mapped path
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Fugaku_RNAseq/Data/Mapped/STAR_mm10_noMultimapFilter/s_GV.WE.Aligned.sortedByCoord.out.bam"

######################################################## READ DATA
# read repeatMasker
rmsk_gr <-
  read_delim(file = rmsk_path, delim = "\t", col_types = cols(start = col_double(), end = col_double())) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

######################################################## MAIN CODE
# ### get count of reads, save summarizedExperiment as RDS
# # load bam file list in memory
# bamfile_chunk <- Rsamtools::BamFile(bam_path, yieldSize = 2000000)
# 
# # register workers for parallel counting
# BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))
# 
# # summarize overlaps
# se <- GenomicAlignments::summarizeOverlaps(features = rmsk_gr,
#                                            reads = bamfile_chunk,
#                                            mode = "Union",
#                                            singleEnd = FALSE,
#                                            ignore.strand = TRUE)
# 
# # save summarizeOverlaps
# saveRDS(se, file = file.path(outpath, "s_GV.WE.repeatMasker.summarizeExperiment.RDS"))

# load summarizeOverlaps
se <- readRDS(file = file.path(outpath, "s_GV.WE.repeatMasker.summarizeExperiment.RDS"))

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

# join with expression counts, summarize per repName
rmsk_repName_count <- 
  rmsk_df %>% 
  dplyr::mutate(GV_count = unlist(assay(se))) %>% 
  dplyr::group_by(repName) %>% 
  dplyr::summarize(GV_count = sum(GV_count) %>% as.numeric(.)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., rmsk_relations, by = "repName") %>% 
  dplyr::left_join(., rmsk_repName_insertions, by = "repName") %>% 
  dplyr::select(repClass, repFamily, repName, insertion_num, GV_count) %>% 
  dplyr::arrange(desc(GV_count)) 

# summarize per repFamily
rmsk_repFamily_count <- 
  rmsk_repName_count %>% 
  dplyr::mutate(repFamily = ifelse(is.na(repFamily), repClass, repFamily)) %>% 
  dplyr::group_by(repFamily) %>% 
  dplyr::summarize(GV_count = sum(GV_count) %>% as.numeric(.), 
                   insertion_num = sum(insertion_num) %>% as.numeric(.), 
                   repClass = str_c(unique(repClass), collapse = "|")) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::select(repClass, repFamily, insertion_num, GV_count) %>% 
  dplyr::arrange(desc(GV_count)) 

# write to separate sheets in xlsx 
write.xlsx(x = rmsk_repName_count %>% as.data.frame(.), 
           file = file.path(outpath, "GV_WE.Fugaku.repeatMasker.counts.xlsx"), 
           sheetName = "repName_sum", 
           row.names = FALSE)
write.xlsx(x = rmsk_repFamily_count %>% as.data.frame(.) , 
           file = file.path(outpath, "GV_WE.Fugaku.repeatMasker.counts.xlsx"), 
           sheetName = "repFamily_sum", 
           append = TRUE, 
           row.names = FALSE)
