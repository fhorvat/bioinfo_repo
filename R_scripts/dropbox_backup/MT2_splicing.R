library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(Rsamtools)
library(BiocParallel)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(tidyr)

rm(list = ls()); gc()
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/splicing")
options(bitmapType = "cairo")

################################################################## functions
findSJwithM2Overlap <- function(x){
  
  # read SJ.out table from STAR
  sjout <- 
    read_delim(sample_df[x, "sjout_path"], delim = "\t", col_names = c("seqnames", "start", "end", "strand", "intron_motif", "annot", "n_uniq", "n_multi", "overhang")) %>%
    mutate(strand = all_strands[as.character(strand)], 
           intron_motif = all_intron_motifs[as.character(intron_motif)], 
           SJ_fullName = paste0(seqnames, ":",
                                start, "-", 
                                end, ":", 
                                strand)) %>%
    filter(annot == 0, intron_motif != "non-canonical") %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
  
  # get coordinates of splice junction start and splice junction end from sj out
  sj_splice_start <- GenomicRanges::resize(sjout, fix = "start", width = 1)
  # sj_splice_start <- c(GenomicRanges::resize(sj_splice_start[strand(sj_splice_start) == "+"], width = 10, fix = "end"), 
  #                      GenomicRanges::resize(sj_splice_start[strand(sj_splice_start) == "-"], width = 10, fix = "start"))
  # sj_splice_start <- GenomicRanges::shift(sj_splice_start, 5)
  
  # overlaps between MT2s and gene splice start
  sj_splice_start_overlapping_MT2 <- sj_splice_start[queryHits(findOverlaps(sj_splice_start, rptmsk_MT2))]
  sj_splice_start_overlapping_MT2$MT2_fullName <- rptmsk_MT2[subjectHits(findOverlaps(sj_splice_start, rptmsk_MT2))]$fullName
  
  # remove those which overlap ENSEMBL/UCSC splice donors
  sj_splice_start_overlapping_MT2 <- sj_splice_start_overlapping_MT2[-queryHits(findOverlaps(sj_splice_start_overlapping_MT2, splice_don_ensembl))]
  sj_splice_start_overlapping_MT2 <- sj_splice_start_overlapping_MT2[-queryHits(findOverlaps(sj_splice_start_overlapping_MT2, splice_don_ucsc))]
  
  # get full splice junctions
  sj_overlapping_MT2 <- sjout[sjout$SJ_fullName %in% sj_splice_start_overlapping_MT2$SJ_fullName]
  sj_overlapping_MT2$MT2_fullName <- sj_splice_start_overlapping_MT2[match(sj_overlapping_MT2$SJ_fullName, sj_splice_start_overlapping_MT2$SJ_fullName)]$MT2_fullName
  sj_overlapping_MT2 <- sj_overlapping_MT2[order(sj_overlapping_MT2$MT2_fullName)]
  
  # final data.frame
  sj_overlapping_MT2_df <- 
    as.data.frame(sj_overlapping_MT2) %>% 
    mutate(MT2_junction_name = paste0(gsub("\\|MT2_Mm", "", MT2_fullName), "|", SJ_fullName), 
           library_size = sample_df[x, "library_size"], 
           RPM_unique = n_uniq / library_size, 
           RPM_multimap = n_multi / library_size) %>% 
    select(MT2_junction_name, RPM_unique, RPM_multimap) %>% 
    arrange(MT2_junction_name)
  colnames(sj_overlapping_MT2_df)[2:3] <- paste0(colnames(sj_overlapping_MT2_df)[2:3], "_", sample_df[x, "sample_name"])
  
  return(sj_overlapping_MT2_df)
}

################################################################## sample data in data.frame
sample_path <- "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2"
sample_df <- 
  data.frame(track_path = list.files(path = sample_path, 
                                     pattern = "*.bam$", 
                                     recursive = T, 
                                     full.names = T), 
             log_path = list.files(path = sample_path, 
                                   pattern = "*Log.final.out", 
                                   recursive = T, 
                                   full.names = T), 
             sjout_path = list.files(path = sample_path, 
                                     pattern = "*SJ.out.tab", 
                                     recursive = T, 
                                     full.names = T), 
             stringsAsFactors = F) %>%
  mutate(sample_name = gsub("^/.*/|\\.bam", "", track_path), 
         library_size = sapply(X = log_path, function(X) as.integer(read.delim(X, header = F, stringsAsFactors = F)[8, 2]) / 10e6)) %>% 
  select(sample_name, track_path, sjout_path, library_size) %>%
  right_join(data.frame(sample_name = c("s_GV.WE", "s_1cell.WE", "s_2cell.WE", "s_2cell.WE_DNAm", "s_4cell.WE")), by = "sample_name")

################################################################## 
# read repeatMasker, get MT2 LTRs, make GRanges

# rptmsk_MT2 <- 
#   read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/coverage/MT2_ORR1A0_solo_LTRs_orderedByFPKMin2cell.txt", delim =  "\t") %>%
#   select(-4) %>%
#   filter(repName == "MT2_Mm") %>%
#   mutate(fullName = paste0(seqnames, ":",
#                            start, "-", 
#                            end, ":", 
#                            strand, "|", 
#                            repName)) %>%
#   arrange(desc(FPKM)) %>%
#   mutate(FPKM_order = dense_rank(desc(FPKM))) %>%
#   makeGRangesFromDataFrame(keep.extra.columns = T)
# 

rptmsk <- read_delim("/common/WORK/fhorvat/reference/mm10/UCSC_repeatMasker_mm10_20161012.txt.gz", delim =  "\t") 

rptmsk_MT2 <- 
  rptmsk %>% 
  filter(element_name == "MT2_Mm") %>%
  mutate(fullName = paste0(seqnames, ":",
                           start, "-", 
                           end, ":", 
                           strand, "|", 
                           element_name)) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

rptmsk <- makeGRangesFromDataFrame(rptmsk, keep.extra.columns = T)

# ################################################################## repeatMasker MT2s FPKM
# # counts over exons of all genes from knownGene table from UCSC
# register(MulticoreParam())
# se_MT2 <- summarizeOverlaps(features = rptmsk_MT2, 
#                             reads = BamFileList(sample_df$track_path, yieldSize = 2000000), 
#                             mode = "Union", 
#                             singleEnd = FALSE, 
#                             ignore.strand = TRUE)
# 
# fpkm_df <- as.data.frame(se_MT2) 
# colnames(fpkm_df) = sample_df$sample_name
# fpkm_df <- 
#   fpkm_df %>%
#   mutate(width = width(rptmsk_MT2), 
#          MT2_fullName = rptmsk_MT2$fullName)
# 
# invisible(lapply(X = sample_df$sample_name,
#                  FUN = function(X) fpkm_df[, X] <<- fpkm_df[, X] / (sample_df[sample_df$sample_name == X, "library_size"] * (fpkm_df$width / 1000))))
# 
# fpkm_df <- 
#   fpkm_df %>%
#   select(-width)
# 
# write.csv(fpkm_df, "fpkm_fugaku_MT2.csv", quote = F, row.names = F)

################################################################## 
# read ENSEMBL splice donors/acceptors
splice_don_ensembl <- 
  read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/splicing/Ensembl_GRCm38.86.20161128_splice_donor.txt", delim =  "\t") %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)
splice_don_ensembl <- c(shift(splice_don_ensembl[strand(splice_don_ensembl) == "+"], 1), 
                        shift(splice_don_ensembl[strand(splice_don_ensembl) == "-"], -1))

splice_acc_ensembl <- 
  read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/splicing/Ensembl_GRCm38.86.20161128_splice_acceptor.txt", delim =  "\t") %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

# read UCSC splice donors/acceptors
splice_don_ucsc <- 
  read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/splicing/UCSC_knownGene_mm10_20161126_splice_donor.txt", delim =  "\t") %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)
splice_don_ucsc <- c(shift(splice_don_ucsc[strand(splice_don_ucsc) == "+"], 1), 
                     shift(splice_don_ucsc[strand(splice_don_ucsc) == "-"], -1))

splice_acc_ucsc <- 
  read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/splicing/UCSC_knownGene_mm10_20161126_splice_acceptor.txt", delim =  "\t") %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

# make non-redundant set of splice acceptors
splice_acc_overlaps <- findOverlaps(splice_acc_ensembl, splice_acc_ucsc)
splice_acc <- c(splice_acc_ensembl[-queryHits(splice_acc_overlaps)], splice_acc_ucsc)
splice_acc <- c(shift(splice_acc[strand(splice_acc) == "+"], -1), 
                shift(splice_acc[strand(splice_acc) == "-"], 1))

################################################################## 
# set strand and intron motifs
all_strands <- setNames(c("*", "+", "-"), as.character(0:2))
all_intron_motifs <- setNames(c("non-canonical", "GT/AG", "CT/AC", "GC/AG", "CT/GC", "AT/AC", "GT/AT"), as.character(0:6))

################################################################## 
# read all splice junctions, get unique set
all_split_junctions <- 
  lapply(X = 1:nrow(sample_df), FUN = function(x){
    sjout <- 
      read_delim(sample_df[x, "sjout_path"], delim = "\t", col_names = c("seqnames", "start", "end", "strand", "intron_motif", "annot", "n_uniq", "n_multi", "overhang")) %>%
      mutate(strand = all_strands[as.character(strand)], 
             intron_motif = all_intron_motifs[as.character(intron_motif)], 
             SJ_fullName = paste0(seqnames, ":",
                                  start, "-", 
                                  end, ":", 
                                  strand)) %>%
      filter(annot == 0, intron_motif != "non-canonical") %>% 
      makeGRangesFromDataFrame(keep.extra.columns = T)
  })

all_split_junctions <- do.call(c, all_split_junctions)
all_split_junctions <- all_split_junctions[!duplicated(all_split_junctions$SJ_fullName)]
mcols(all_split_junctions) <- mcols(all_split_junctions)[c(1, 6)]

# get all split junction ends (for downstream splicing annotation)
all_split_junctions_ends <- GenomicRanges::resize(all_split_junctions, fix = "end", width = 1)
# all_split_junctions_ends <- c(GenomicRanges::resize(all_split_junctions_ends[strand(all_split_junctions_ends) == "+"], width = 10, fix = "start"),
#                               GenomicRanges::resize(all_split_junctions_ends[strand(all_split_junctions_ends) == "-"], width = 10, fix = "end"))
# sj_splice_end <- GenomicRanges::shift(all_split_junctions_ends, -5)

# find overlaps of splice ends with annotated exons
all_split_ends_exons_overlaps <- findOverlaps(all_split_junctions_ends, splice_acc)
all_split_ends_exons_overlaps_exons <- splice_acc[subjectHits(all_split_ends_exons_overlaps)]
all_split_ends_exons_overlaps_sj <- 
  as.data.frame(all_split_junctions_ends[queryHits(all_split_ends_exons_overlaps)]) %>% 
  mutate(downstream_splice_end = paste0(all_split_ends_exons_overlaps_exons$transcript_id, ":",
                                        all_split_ends_exons_overlaps_exons$ex.num)) %>% 
  select(SJ_fullName, downstream_splice_end) 

# find overlaps of splice ends with repeatMasker
all_split_ends_rptmsk_overlaps <- findOverlaps(all_split_junctions_ends, rptmsk, ignore.strand = T)
all_split_ends_rptmsk_overlaps_rptmsk <- rptmsk[subjectHits(all_split_ends_rptmsk_overlaps)]
all_split_ends_rptmsk_overlaps_sj <- 
  as.data.frame(all_split_junctions_ends[queryHits(all_split_ends_rptmsk_overlaps)]) %>% 
  mutate(downstream_splice_end = paste0(all_split_ends_rptmsk_overlaps_rptmsk$element_name, ":", strand(all_split_ends_rptmsk_overlaps_rptmsk))) %>% 
  select(SJ_fullName, downstream_splice_end) 

# merge downstream splice site data with splice junctions
all_split_junctions <- 
  as.data.frame(all_split_junctions) %>%
  left_join(all_split_ends_exons_overlaps_sj, by = "SJ_fullName") %>% 
  left_join(all_split_ends_rptmsk_overlaps_sj, by = "SJ_fullName") %>% 
  mutate(downstream_splice = gsub("^NA\\||\\|NA$", "", paste0(downstream_splice_end.x, "|", downstream_splice_end.y))) %>% 
  select(-c(downstream_splice_end.x, downstream_splice_end.y))

################################################################## 
# get list of data frames with MT2s overlaping split junctions
sj_overlapping_MT2_list <- lapply(1:nrow(sample_df), findSJwithM2Overlap)
names(sj_overlapping_MT2_list) <- sample_df$sample_name

# get MT2 elements FPKMs
MT2_fpkm <-
  read_csv("fpkm_fugaku_MT2.csv") %>% 
  select(MT2_fullName, MT2_fpkm_GV = s_GV.WE, MT2_fpkm_1C = s_1cell.WE, 
         MT2_fpkm_2C = s_2cell.WE, MT2_fpkm_2C_DNAm = s_2cell.WE_DNAm, MT2_fpkm_4C = s_4cell.WE)
  
# merge all samples, add data, write as .csv 
sj_overlapping_MT2_all <-
  Reduce(function(...) merge(..., by = "MT2_junction_name", all.x = TRUE, all.y = TRUE), sj_overlapping_MT2_list) %>%
  separate(MT2_junction_name, c("MT2_fullName", "SJ_fullName"), "\\|") %>% 
  mutate(MT2_fullName = paste0(MT2_fullName, "|MT2_Mm")) 
colnames(sj_overlapping_MT2_all) <- gsub("s_|\\.WE", "", colnames(sj_overlapping_MT2_all))
colnames(sj_overlapping_MT2_all) <- gsub("cell", "C", colnames(sj_overlapping_MT2_all))

sj_overlapping_MT2_all <- 
  sj_overlapping_MT2_all %>% 
  left_join(all_split_junctions[, c("SJ_fullName", "intron_motif", "downstream_splice")], by = "SJ_fullName") %>%
  left_join(MT2_fpkm, by = "MT2_fullName") %>% 
  select(MT2_id = MT2_fullName, SJ_id = SJ_fullName, intron_motif, downstream_splice,
         matches("GV"), matches("1C"), matches("2C$"), matches("2C_DNAm"), matches("4C"))

write.csv(sj_overlapping_MT2_all, "MT2_spliced_merged.csv", quote = F, row.names = F)
