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
library(ggthemes)

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

# library size path
library_size_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/library_size", str_c("library_size.counts.", experiment, ".csv"))

# read library size table
library_size_df <- 
  readr::read_csv(library_size_path) %>% 
  tidyr::gather(mismatch, library_size, -sample_id) %>% 
  tidyr::separate(mismatch, c("mismatch", "read_group"), sep = "\\.") %>% 
  tidyr::unite(sample_id, sample_id, mismatch, sep = ".") %>% 
  tidyr::spread(read_group, library_size) %>% 
  dplyr::select(sample_id, library_size = `21to23nt_reads`)


#####
# get counts of MosIR from all bam files in one data.frame, sum 0 mismatches with 1 and 2 mismatches counts
counts_df <- 
  purrr::map(bam_paths, countIR) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::left_join(., library_size_df, by = "sample_id") %>%
  dplyr::mutate(library_size = round((library_size / 1e6), 4),
                alignment_count = alignment_count / library_size) %>% 
  dplyr::select(-library_size) %>% 
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

# #### write
# # write to separate sheets in xlsx 
# write.xlsx(x = lengths_sum %>% 
#              dplyr::select_at(vars("sample_id", starts_with("mis_0"))) %>% 
#              magrittr::set_colnames(., str_replace(colnames(.), "mis_[0-2]\\.", "l_")) %>% 
#              as.data.frame(.), 
#            file = file.path(outpath, str_c("IR.", experiment, ".read_length.RPM.xlsx")), 
#            sheetName = "allowed_0_mismatches", 
#            row.names = FALSE)
# 
# write.xlsx(x = lengths_sum %>% 
#              dplyr::select_at(vars("sample_id", starts_with("mis_1"))) %>% 
#              magrittr::set_colnames(., str_replace(colnames(.), "mis_[0-2]\\.", "l_")) %>% 
#              as.data.frame(.), 
#            file = file.path(outpath, str_c("IR.", experiment, ".read_length.RPM.xlsx")), 
#            sheetName = "allowed_1_mismatches", 
#            append = TRUE,
#            row.names = FALSE)
# 
# write.xlsx(x = lengths_sum %>% 
#              dplyr::select_at(vars("sample_id", starts_with("mis_2"))) %>% 
#              magrittr::set_colnames(., str_replace(colnames(.), "mis_[0-2]\\.", "l_")) %>% 
#              as.data.frame(.), 
#            file = file.path(outpath, str_c("IR.", experiment, ".read_length.RPM.xlsx")), 
#            sheetName = "allowed_2_mismatches",
#            append = TRUE, 
#            row.names = FALSE)
# 
# # summary of all reads
# plasmid_sum <- 
#   counts_df %>%
#   dplyr::group_by(sample_id, mismatch) %>% 
#   dplyr::summarize(alignment_count = sum(alignment_count)) %>% 
#   dplyr::ungroup(.) %>% 
#   tidyr::spread(mismatch, alignment_count) %>% 
#   dplyr::mutate_at(vars(starts_with("mis")), funs(replace(., is.na(.), 0))) %>%
#   dplyr::mutate(mis_1 = mis_1 + mis_0,
#                 mis_2 = mis_2 + mis_1) %>%
#   dplyr::arrange(sample_id)
# 
# # write to separate sheets in xlsx 
# write.xlsx(x = plasmid_sum %>% 
#              as.data.frame(.), 
#            file = file.path(outpath, str_c("IR.", experiment, ".read_length.RPM.xlsx")), 
#            sheetName = "all_reads_summary",
#            append = TRUE,
#            row.names = FALSE)

### plot
plot_df <- 
  lengths_sum %>% 
  dplyr::select_at(vars("sample_id", starts_with("mis_0"))) %>% 
  tidyr::gather(., key = alignment_length, value = RPM, -sample_id) %>% 
  dplyr::mutate(alignment_length = str_remove(alignment_length, "mis_0."), 
                sample_id = str_remove_all(sample_id, "s_|_r[1,2]{1}.SE")) %>% 
  dplyr::filter(alignment_length >= 19, 
                alignment_length <= 24, 
                !str_detect(sample_id, "Elav")) %>% 
  dplyr::group_by(sample_id, alignment_length) %>% 
  dplyr::summarise(avg_RPM = round(mean(RPM)), 
                   sd = sd(RPM)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(sample_id = factor(sample_id, levels = c("pCag_EGFP_mosIR", "pU6_MosIR", "pCag_EGFP_Lin28IR", "pU6_Lin28IR")))

# barplot with error bars
ggplot(data = plot_df, mapping = aes(alignment_length, avg_RPM, fill = sample_id)) +
  geom_bar(stat = "identity", position = position_dodge2(width = 1)) +
  geom_errorbar(aes(ymin = avg_RPM - sd, ymax = avg_RPM + sd),
                width = 0.2, position = position_dodge(0.9), color = "gray30") + 
  scale_y_continuous(breaks = seq(100, 1400, 100)) +
  scale_fill_manual(values = c("#000000", "#7f7f7f", "#4472c4", "#9dc3e6")) + 
  guides(fill = guide_legend(nrow = 1), 
         color = F) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 15, vjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(), 
        legend.text = element_text(size = 15),
        panel.grid.major.x = element_blank(), 
        legend.position = "bottom") +
  ggsave(file.path(outpath, "test2.png"), width = 15, height = 10)


