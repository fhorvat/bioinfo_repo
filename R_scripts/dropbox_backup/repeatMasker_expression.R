library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(magrittr)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(BiocParallel)


setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/expression")

# sample data.frame
sample_df <- 
  data.frame(track_path = list.files(path = "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2", 
                                     pattern = "*.bam$", 
                                     recursive = T, 
                                     full.names = T), 
             log_path = list.files(path = "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2", 
                                   pattern = "*Log.final.out", 
                                   recursive = T, 
                                   full.names = T), 
             sjout_path = list.files(path = "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2", 
                                     pattern = "*SJ.out.tab", 
                                     recursive = T, 
                                     full.names = T), 
             stringsAsFactors = F) %>%
  mutate(sample_name = gsub("^/.*/|\\.bam", "", track_path), 
         library_size = sapply(X = log_path, function(X) as.integer(read.delim(X, header = F, stringsAsFactors = F)[8, 2]) / 10e6)) %>% 
  select(sample_name, track_path, sjout_path, library_size) 

# repeat masker
rptmsk <- read_delim("/common/WORK/fhorvat/reference/mm10/UCSC_repeatMasker_mm10_20161012.txt.gz", delim =  "\t")

rptmsk_LTR <-
  rptmsk %>% 
  filter(element_class == "LTR") %>%
  mutate(fullName = paste0(seqnames, ":",
                           start, "-", 
                           end, ":", 
                           strand, "|", 
                           element_name)) %>% 
  dplyr::select(-element_class) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

# counts over repeatMasker LTRs
register(MulticoreParam())
se_LTR <- summarizeOverlaps(features = rptmsk_LTR,
                            reads = BamFileList(sample_df$track_path, yieldSize = 2000000),
                            mode = "Union",
                            singleEnd = FALSE,
                            ignore.strand = TRUE)

# FPKM
fpkm_df <- as.data.frame(assay(se_LTR))
colnames(fpkm_df) <- sample_df$sample_name
fpkm_df %<>%
  mutate(width = width(rptmsk_LTR),
         LTR_fullName = rptmsk_LTR$fullName)

invisible(lapply(X = sample_df$sample_name,
                 FUN = function(X) fpkm_df[, X] <<- fpkm_df[, X] / (sample_df[sample_df$sample_name == X, "library_size"] * (fpkm_df$width / 1000))))

fpkm_df %<>%
  dplyr::select(-width) %>% 
  write_csv("repeatMasker_LTR_FPKM_Fugaku.csv")
