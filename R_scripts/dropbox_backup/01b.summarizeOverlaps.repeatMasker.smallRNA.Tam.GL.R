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
# bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Fugaku_RNAseq/Data/Mapped/STAR_mm10_noMultimapFilter/s_GV.WE.Aligned.sortedByCoord.out.bam"
bam_path <- 
  list.files(path = c("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/GarciaLopez_2015_RNA_GSE59254/Data/Mapped/STAR_mm10_new", 
                      "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/Tam_2008_Nature_GSE10364/Data/Mapped/STAR_mm10_new"),
             pattern = ".*21to23nt.bam$|.*24to30nt.bam$", 
             full.names = T, recursive = T) %>% 
  .[!str_detect(., "sperm|PGC")]

# GV count path
GV_path <- file.path(outpath, "GV_WE.Fugaku.repeatMasker.counts.xlsx")

######################################################## READ DATA
# read repeatMasker
rmsk_gr <-
  read_delim(file = rmsk_path, delim = "\t", col_types = cols(start = col_double(), end = col_double())) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

# read GV counts - repName
GV_repName_count <- 
  xlsx::read.xlsx(file = GV_path, sheetName = "repName_sum") %>% 
  dplyr::select(repName, GV_count) %>% 
  dplyr::mutate_at(vars(starts_with("rep")), funs(as.character(.)))

# read GV counts - repFamily
GV_repFamily_count <- 
  xlsx::read.xlsx(file = GV_path, sheetName = "repFamily_sum") %>% 
  dplyr::select(repFamily, GV_count) %>% 
  dplyr::mutate_at(vars(starts_with("rep")), funs(as.character(.)))

######################################################## MAIN CODE
# ### get count of reads, save summarizedExperiment as RDS
# # load bam file list in memory
# bamfile_chunk <- Rsamtools::BamFileList(bam_path, yieldSize = 2000000)
# 
# # register workers for parallel counting
# BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))
# 
# # summarize overlaps
# se <- GenomicAlignments::summarizeOverlaps(features = rmsk_gr,
#                                            reads = bamfile_chunk,
#                                            mode = "Union",
#                                            singleEnd = TRUE,
#                                            ignore.strand = TRUE)
# 
# # save summarizeOverlaps
# saveRDS(se, file = file.path(outpath, "smallRNA.Tam.GL.repeatMasker.summarizeExperiment.RDS"))

# load summarizeOverlaps
se <- readRDS(file = file.path(outpath, "smallRNA.Tam.GL.repeatMasker.summarizeExperiment.RDS"))

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
  dplyr::select_at(.vars = vars(matches("s_MII.SE.*|s_oocyte_19to30.*"))) %>% 
  magrittr::set_colnames(., str_remove_all(colnames(.), ".bam|.SE") %>% 
                           str_replace(., "s_MII", "GarciaLopez.MII") %>% 
                           str_replace(., "s_oocyte_", "Tam.oocyte."))

# join with expression counts, summarize per repName
rmsk_repName_count <- 
  rmsk_df %>% 
  dplyr::bind_cols(., se_df) %>% 
  dplyr::group_by(repName) %>% 
  dplyr::summarize_at(.vars = vars(matches("GarciaLopez.*|Tam.*")), 
                      .funs = funs(sum(.) %>% as.numeric(.))) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., rmsk_relations, by = "repName") %>% 
  dplyr::left_join(., rmsk_repName_insertions, by = "repName") %>% 
  dplyr::left_join(., GV_repName_count, by = "repName") %>% 
  dplyr::select(repClass, repFamily, repName, insertion_num, GV_count, everything()) %>% 
  dplyr::arrange(desc(GV_count))

# summarize per repFamily
rmsk_repFamily_count <- 
  rmsk_repName_count %>% 
  dplyr::mutate(repFamily = ifelse(is.na(repFamily), repClass, repFamily)) %>% 
  dplyr::group_by(repFamily)

rmsk_repFamily_count <-
  left_join(x = rmsk_repFamily_count %>% dplyr::summarize_at(vars(matches("GarciaLopez.*|Tam.*|insertion_num")), funs(sum(.) %>% as.numeric(.))),
            y = rmsk_repFamily_count %>% dplyr::summarise(repClass = str_c(unique(repClass), collapse = "|")), 
            by = "repFamily") %>% 
  dplyr::left_join(., GV_repFamily_count, by = "repFamily") %>% 
  dplyr::select(repClass, repFamily, insertion_num, GV_count, everything()) %>% 
  dplyr::arrange(desc(GV_count))

# write to separate sheets in xlsx 
write.xlsx(x = rmsk_repName_count %>% as.data.frame(.), 
           file = file.path(outpath, "smallRNA.Tam.GL.repeatMasker.counts.xlsx"), 
           sheetName = "repName_sum", 
           row.names = FALSE)
write.xlsx(x = rmsk_repFamily_count %>% as.data.frame(.) , 
           file = file.path(outpath, "smallRNA.Tam.GL.repeatMasker.counts.xlsx"), 
           sheetName = "repFamily_sum", 
           append = TRUE, 
           row.names = FALSE)
