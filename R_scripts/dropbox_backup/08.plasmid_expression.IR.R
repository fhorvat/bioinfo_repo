### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/IR_expression")

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

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
### get counts of read alignments with different mismatch cut-off over MosIR in plasmid 
countIR <- function(bam_path){
  
  # get bam name
  bam_name <- basename(bam_path) %>% str_remove(., ".bam")
  
  # get IR name
  ir_name <- str_remove_all(bam_name, "^s_|_r[1,2].SE.mis_[0-2]") %>% str_replace(., "pCag_EGFP", "pCag-EGFP")
  
  # subset IR GenomicRanges
  ir_filt <- ir_gr[seqnames(ir_gr) == ir_name]
  
  # read bam
  bam_gr <- 
    GenomicAlignments::readGAlignmentsList(file = bam_path, 
                                           use.names = TRUE,
                                           param = ScanBamParam(which = ir_filt, 
                                                                flag = scanBamFlag(isMinusStrand = NA))) %>% 
    unlist(.)
  
  # return tibble
  count_table <- 
    tibble(read_id = names(bam_gr), 
           alignment_length = width(bam_gr)) %>% 
    dplyr::group_by(read_id) %>% 
    dplyr::summarise(alignment_length = max(alignment_length)) %>% 
    dplyr::group_by(alignment_length) %>% 
    dplyr::summarise(alignment_count = n()) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::right_join(tibble(alignment_length = rep(12:50)), by = "alignment_length") %>% 
    dplyr::mutate(alignment_count = replace(alignment_count, is.na(alignment_count), 0)) %>% 
    dplyr::mutate(sample_id = bam_name)
  
  # return
  return(count_table)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# get path of inverted repeat coordinates from mRNA
ir_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/mm10.all_plasmids/IR_mRNAs.plasmid_coordinates/IR_mRNA.plasmid_coordinates.clean.csv"

######################################################## READ DATA
# read coordinates of inverted repeats from mRNA
ir_gr <- readr::read_csv(file = ir_path)

######################################################## MAIN CODE
# prepare mRNA IR coordinates on plasmids
ir_gr %<>% GenomicRanges::GRanges(.)
mcols(ir_gr) <- mcols(ir_gr)[, c("arm")]
names(mcols(ir_gr)) <- "gene_id"
mcols(ir_gr)$gene_biotype <- as.character(seqnames(ir_gr))

### get IR counts for all samples in each experiment 
# experiment paths
experiment_path <- 
  "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/NIH3T3_transfected.2018/Data/Mapped/STAR_mm10.plasmids" %>% 
  magrittr::set_names(., str_remove(., "Data.*$") %>% basename(.))

# list of experiments
experiment_list <- names(experiment_path)

# loop through experiments
experiment <- experiment_list[1]

# get mapped path 
mapped_path <- experiment_path[experiment]

# list bams 
bam_paths <- list.files(path = mapped_path, pattern = ".*mis_.*bam$", full.names = T)

# get counts of MosIR from all bam files in one data.frame, sum 0 mismatches with 1 and 2 mismatches counts
counts_df <- 
  purrr::map(bam_paths, countIR) %>% 
  dplyr::bind_rows(.) %>% 
  tidyr::separate(sample_id, c("sample_id", "mismatch"), sep = ".mis_") %>% 
  dplyr::mutate(mismatch = str_c("mis_", mismatch))

# summary of reads of different lengths on one plasmid
lengths_sum <- 
  counts_df %>% 
  tidyr::spread(mismatch, alignment_count) %>% 
  dplyr::mutate_at(vars(starts_with("mis")), funs(replace(., is.na(.), 0))) %>% 
  dplyr::mutate(mis_1 = mis_1 + mis_0, 
                mis_2 = mis_2 + mis_1) %>% 
  dplyr::arrange(sample_id, alignment_length) %>% 
  dplyr::filter(alignment_length >= 18, alignment_length <= 30) %>% 
  tidyr::gather(mismatch_allowed, alignment_count, mis_0:mis_2) %>% 
  tidyr::unite(temp, mismatch_allowed, alignment_length, sep = ".") %>%
  tidyr::spread(temp, alignment_count) %>%
  dplyr::mutate_at(vars(starts_with("mis_")), funs(replace(., is.na(.), 0)))

# write to separate sheets in xlsx 
write.xlsx(x = lengths_sum %>% 
             dplyr::select_at(vars("sample_id", starts_with("mis_0"))) %>% 
             magrittr::set_colnames(., str_replace(colnames(.), "mis_[0-2]\\.", "l_")) %>% 
             as.data.frame(.), 
           file = file.path(outpath, str_c("IR.", experiment, ".read_length.counts.xlsx")), 
           sheetName = "allowed_0_mismatches", 
           row.names = FALSE)

write.xlsx(x = lengths_sum %>% 
             dplyr::select_at(vars("sample_id", starts_with("mis_1"))) %>% 
             magrittr::set_colnames(., str_replace(colnames(.), "mis_[0-2]\\.", "l_")) %>% 
             as.data.frame(.), 
           file = file.path(outpath, str_c("IR.", experiment, ".read_length.counts.xlsx")), 
           sheetName = "allowed_1_mismatches", 
           append = TRUE,
           row.names = FALSE)

write.xlsx(x = lengths_sum %>% 
             dplyr::select_at(vars("sample_id", starts_with("mis_2"))) %>% 
             magrittr::set_colnames(., str_replace(colnames(.), "mis_[0-2]\\.", "l_")) %>% 
             as.data.frame(.), 
           file = file.path(outpath, str_c("IR.", experiment, ".read_length.counts.xlsx")), 
           sheetName = "allowed_2_mismatches",
           append = TRUE, 
           row.names = FALSE)

# summary of all reads
plasmid_sum <- 
  counts_df %>%
  dplyr::group_by(sample_id, mismatch) %>% 
  dplyr::summarize(alignment_count = sum(alignment_count)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::spread(mismatch, alignment_count) %>% 
  dplyr::mutate_at(vars(starts_with("mis")), funs(replace(., is.na(.), 0))) %>%
  dplyr::mutate(mis_1 = mis_1 + mis_0,
                mis_2 = mis_2 + mis_1) %>%
  dplyr::arrange(sample_id)

# write to separate sheets in xlsx 
write.xlsx(x = plasmid_sum %>% 
             as.data.frame(.), 
           file = file.path(outpath, str_c("IR.", experiment, ".read_length.counts.xlsx")), 
           sheetName = "all_reads_summary",
           append = TRUE,
           row.names = FALSE)


