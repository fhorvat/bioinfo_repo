### INFO: counts of multimap reads 
### DATE: 31. 08. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()
# options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Fugaku_RNAseq/Analysis/multimap_reads/rmsk_reads_classification")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tibble)
library(data.table)
library(doMC)

library(Rsamtools)
library(Biostrings)
library(GenomicRanges)
library(GenomicAlignments)
library(BiocParallel)

######################################################## PATH VARIABLES
outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/Fugaku_RNAseq/Analysis/multimap_reads/rmsk_reads_classification/results"
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Fugaku_RNAseq/Data/Mapped/STAR_mm10_noMultimapFilter"
bam_list <- list.files(path = bam_path, pattern = "*bam$", recursive = F, full.names = T)

lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
repeatmasker_path <- "/common/WORK/fhorvat/reference/mouse/mm10/UCSC/UCSC_repeatMasker_withLTRsubclass.txt.gz"

file <- bam_list[1]
file_name <- str_replace_all(file, "\\/.*\\/|_Aligned.sortedByCoord.out.bam", "")

######################################################## SOURCE FILES
source(file.path(lib_path, "GffToGRanges.R"))
source(file.path(lib_path, "headt.R"))
source("/common/WORK/fhorvat/code_library/R_scripts/vfranke/BamWorkers.R")

######################################################## FUNCTIONS

######################################################## READ DATA
# repeatMaskerVIZ - retrotransposons only
rptmsk_gr <-
  read_delim(file = repeatmasker_path, delim = "\t") %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

# get LTR family data
ltr_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/LTR_genomic_insertion_freq/LTR_families_data.csv"
LTR_data <- readr::read_csv(file = ltr_path)

######################################################## MAIN CODE
### read multimappers (NH > 1) from .bam
# reads in the .bam
bam <- suppressWarnings(readGAlignmentPairs(file, param = ScanBamParam(tagFilter = list(NH = c(2:999))), use.names = TRUE))
bams <- readGAlignments(file, param = ScanBamParam(tagFilter = list(NH = c(2:999)), tag = "NH"), use.names = TRUE)

# get coordinates of each read alignment
gbam <- grglist(bam)

# get read ID's and number of alignment position (NH tag) 
tab <- 
  tibble(read_id = names(bams), NH = values(bams)$NH) %>% 
  dplyr::distinct(.) %>% 
  dplyr::filter(read_id %in% names(bam))

rm(bam, bams); gc()


### get reads and align class (LTR > LINE > SINE > other > non-repeat)
# get reads which align to repeats
reads_all <- 
  findOverlaps(gbam, rptmsk_gr, ignore.strand = TRUE) %>% 
  as.matrix(.) %>% 
  tibble::as.tibble(.) %>% 
  dplyr::mutate(read_id = names(gbam)[queryHits], 
                repClass = rptmsk_gr$repClass[subjectHits],
                LTR_subclass = rptmsk_gr$LTR_subclass[subjectHits]) %>% 
  dplyr::select(-(1:2)) 

# get read ID's of all reads which align to LTRs of interest
reads_LTR <- 
  reads_all %>% 
  dplyr::filter(!is.na(LTR_subclass)) %>% 
  dplyr::distinct(.) %>% 
  dplyr::group_by(read_id) %>% 
  dplyr::summarise(LTR_subclass = stringr::str_c(LTR_subclass, collapse = "|")) %>% 
  dplyr::mutate(LTR_subclass = replace(LTR_subclass, stringr::str_detect(LTR_subclass, "\\|"), "LTR_mixed")) %>% 
  dplyr::rename(final_class = LTR_subclass)

# get reads which map to other repeats
reads_all %<>% 
  dplyr::filter(!(read_id %in% reads_LTR$read_id)) %>% 
  dplyr::select(-LTR_subclass) %>% 
  dplyr::distinct(.) %>% 
  dplyr::mutate(final_class = NA)

# sequentialy add class to reads aligning to different repeat types
for(repeat_filter in c("LTR", "LINE", "SINE")){
  
  # get ID's of read 
  filter_id <- 
    reads_all %>% 
    dplyr::filter(repClass == repeat_filter & is.na(final_class)) %$% 
    read_id %>% 
    unique(.)
  
  # add class to the read
  reads_all %<>%
    dplyr::mutate(final_class = ifelse(read_id %in% filter_id, repeat_filter, final_class))
  
}

# get distinct read ID's-class combinations, join with LTR reads, join with reads which don't align to repeatMasker
reads_all %<>%
  dplyr::mutate(final_class = replace(final_class, is.na(final_class), "repeat_other"), 
                final_class = replace(final_class, final_class == "LTR", "LTR_other")) %>% 
  dplyr::select(-repClass) %>% 
  dplyr::distinct(.) %>% 
  rbind(reads_LTR) %>% 
  rbind(tibble(read_id = unique(names(gbam)[!(names(gbam) %in% .$read_id)]), 
               final_class = "not_repeat")) %>% 
  dplyr::left_join(tab, by = "read_id") 


### plot repeat classes
# data frame for plot
reads_class_plot <- 
  reads_all %>% 
  dplyr::mutate(final_class = replace(final_class, final_class %in% unique(reads_LTR$final_class), "LTR_interest"), 
                final_class = factor(final_class, levels = c("LTR_interest", "LTR_other", "LINE", "SINE", "repeat_other", "not_repeat"))) %>% 
  dplyr::group_by(NH) %>% 
  dplyr::count(final_class) %>% 
  dplyr::arrange(desc(n))

# histogram
ggplot(reads_class_plot) +
  geom_rect(aes(xmin = NH, ymin = 0, xmax = NH + 1, ymax = n, fill = final_class), color = "black") +
  scale_x_continuous(limits = c(0, 100)) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_blank()) +
  facet_grid(final_class ~ ., scales = "free") +
  ggsave(filename = file.path(outpath, stringr::str_c(file_name, "_multimap_histogram_repeats.pdf")), 
         width = 30, height = 20, units = "cm")

# zoomed histogram
ggplot(reads_class_plot) +
  geom_rect(aes(xmin = NH, ymin = 0, xmax = NH + 1, ymax = n, fill = final_class), color = "black") +
  scale_x_continuous(limits = c(0, 100)) +
  theme_bw() + 
  coord_cartesian(ylim = c(0, 1000)) +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_blank()) +
  facet_grid(final_class ~ ., scales = "free") +
  ggsave(filename = file.path(outpath, stringr::str_c(file_name, "_multimap_histogram_repeats_zoom.pdf")), 
         width = 30, height = 20, units = "cm")


### plot LTR classes
# order LTR's by multimapper counts
LTR_topN <-
  reads_LTR %>% 
  dplyr::count(final_class) %>%
  dplyr::arrange(desc(n)) %>% 
  dplyr::top_n(n = 15, wt = n) %$% 
  final_class

# data frame for plot
LTR_class_plot <- 
  reads_LTR %>% 
  dplyr::filter(final_class %in% LTR_topN) %>% 
  dplyr::left_join(tab, by = "read_id") %>%
  dplyr::group_by(NH) %>% 
  dplyr::count(final_class) %>% 
  dplyr::mutate(final_class = factor(final_class, levels = LTR_topN))

# histogram plot
ggplot(LTR_class_plot) +
  geom_rect(aes(xmin = NH, ymin = 0, xmax = NH + 1, ymax = n, fill = final_class), color = "black") +
  scale_x_continuous(limits = c(0, 100)) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_blank()) +
  facet_grid(final_class ~ ., scales = "free") +
  ggsave(filename = file.path(outpath, stringr::str_c(file_name, "_multimap_histogram_LTRs_top.pdf")), 
         width = 30, height = 50, units = "cm")

### all LTR's
# order LTR's by multimapper counts
LTR_topN <-
  reads_LTR %>% 
  dplyr::count(final_class) %>%
  dplyr::arrange(desc(n)) %$% 
  final_class

# data frame for plot
LTR_class_plot <- 
  reads_LTR %>% 
  dplyr::left_join(tab, by = "read_id") %>%
  dplyr::group_by(NH) %>% 
  dplyr::count(final_class) %>% 
  dplyr::mutate(final_class = factor(final_class, levels = LTR_topN))

# histogram plot
ggplot(LTR_class_plot) +
  geom_rect(aes(xmin = NH, ymin = 0, xmax = NH + 1, ymax = n, fill = final_class), color = "black") +
  scale_x_continuous(limits = c(0, 100)) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_blank()) +
  facet_grid(final_class ~ ., scales = "free") +
  ggsave(filename = file.path(outpath, stringr::str_c(file_name, "_multimap_histogram_LTRs_all.pdf")), 
         width = 30, height = 100, units = "cm")

# ### density plot
# # data frame
# reads_class_plot <- 
#   reads_all %>% 
#   dplyr::mutate(final_class = replace(final_class, final_class %in% unique(reads_LTR$final_class), "LTR_interest"), 
#                 final_class = factor(final_class, levels = c("LTR_interest", "LTR_other", "LINE", "SINE", "repeat_other", "not_repeat")))
# 
# # plot
# ggplot(reads_class_plot, aes(NH, ..count.., color = final_class)) +
#   stat_density(geom = "line") +
#   scale_x_continuous(limits = c(0, 100)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks.x = element_blank()) +
#   facet_grid(final_class ~ ., scales = "free") +
#   ggsave(filename = file.path(outpath, stringr::str_c(file_name, "_multimap_density.pdf")),
#          width = 30, height = 20, units = "cm")