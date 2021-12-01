library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(magrittr)
library(readr)
library(stringr)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(BiocParallel)

# options(bitmapType = "cairo")
# rm(list = ls()[!grepl("bam_|coverage_", ls())])

################################################################################## functions

################################################################################## reading data
# setting working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/cow_expression")

# list of chosen cow LTRs
chosen_LTRs <- c("ERV1-1-LTR_BT", "ERV1-2-LTR_BT", "LTR32_BT", "MLT1D",
                 "MLT1A0", "BTLTR1", "ERV1-1B-LTR_BT", "MLT1B", "MLT1C",
                 "MER41_BT", "LTR75_BT", "LTR31B_BT", "MER21C", "LTR16A",
                 "MLT1K", "MER87B_BT", "MLT1A", "MLT1E3", "MLT1J", "ERV1-3-LTR_BT",
                 "ERV54-EC_LTR", "LTR13B_BT", "LTR33", "LTR50", "LTR79",
                 "MLT1F", "MLT1G3", "MLT1H", "MLT1I", "MLT1L",
                 "ALTR2C_BT", "ERV1-2C-LTR_BT", "ERV2-1C-LTR_BT")

# making TxDb object from refseq gtf
refseq_genes <- 
  makeTxDbFromGFF("/common/WORK/fhorvat/reference/cow/bosTau7/160815.bosTau7.UCSC.Refseq.gtf") %>% 
  genes(.)

# repeatMasker table
repeatMasker <- 
  read_delim("/common/WORK/fhorvat/reference/cow/bosTau7/UCSC_repeatMasker_bosTau7_20170206_all_fields.txt.gz", delim = "\t") %>% 
  dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, repName, repClass, repFamily)

# take only chosen LTR
repeatMasker_filtered <- 
  repeatMasker %>% 
  dplyr::filter(repName %in% chosen_LTRs) %>% 
  dplyr::select(-c(repClass, repFamily)) %>% 
  mutate(fullName = paste0(seqnames, ":", 
                           start, "-", 
                           end, ":",
                           strand, "|", 
                           repName)) %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T)

################################################################################## .bam files
# .bam files paths
filenames <- file.path(c("/common/DB/vfranke/Base/AccessoryData/Graf_2014_PNAS/Mapped/STAR_bosTau7_Merged/s_GV/s_GV.bam",
                         "/common/DB/vfranke/Base/AccessoryData/Graf_2014_PNAS/Mapped/STAR_bosTau7_Merged/s_MII/s_MII.bam",
                         "/common/DB/vfranke/Base/AccessoryData/Graf_2014_PNAS/Mapped/STAR_bosTau7_Merged/s_4c/s_4c.bam",
                         "/common/DB/vfranke/Base/AccessoryData/Graf_2014_PNAS/Mapped/STAR_bosTau7_Merged/s_8c/s_8c.bam",
                         "/common/DB/vfranke/Base/AccessoryData/Graf_2014_PNAS/Mapped/STAR_bosTau7_Merged/s_16c/s_16c.bam",
                         "/common/DB/vfranke/Base/AccessoryData/Graf_2014_PNAS/Mapped/STAR_bosTau7_Merged/s_Blast/s_Blast.bam"))

## get library size in million of reads
library_size_df <- 
  read_delim(file = "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/cow_expression/Graf_2014_library_size.txt", delim = "\t") 

######################################################################### getting count of reads over repeatMasker filtered table
# counting overlaps
bamfiles <- BamFileList(filenames, yieldSize = 2000000)
register(MulticoreParam(workers = 8))
se <- summarizeOverlaps(features = repeatMasker_filtered,
                        reads = bamfiles,
                        mode = "Union",
                        singleEnd = TRUE,
                        ignore.strand = TRUE)

# RPM from counts per selected repNames
rpm_df <-
  assay(se) %>% 
  as.data.frame() %>% 
  set_colnames(str_replace(colnames(.), ".bam", "")) %>% 
  mutate(repName = repeatMasker_filtered$repName) %>% 
  dplyr::select(repName, 1:6) %>% 
  group_by(repName) %>% 
  summarise_all(.funs = sum) %>% 
  tidyr::gather(key = sample, value = rpm, -repName) %>% 
  left_join(library_size_df, by = "sample") %>% 
  mutate(rpm = rpm / library_size, 
         rpm = log10(rpm)) %>% 
  dplyr::select(-library_size) %>% 
  mutate(sample = factor(sample, levels = c("s_GV", "s_MII", "s_4c", "s_8c", "s_16c", "s_Blast")), 
         repName = factor(repName, levels = rev(unique(repName)))) 
  
# plot as heatmap
ggplot(data = rpm_df, aes(x = sample, y = repName)) + 
  geom_tile(aes(fill = rpm), colour = "white") + 
  scale_fill_gradient(low = "white", high = "black") +
  guides(fill = guide_colorbar(ticks = FALSE, 
                               title = expression(paste(log[10], " RPM")), 
                               title.position = "left", 
                               label.position = "left")) + 
  scale_x_discrete(position = "top", 
                   labels = c("GV", "MII", "4C", "8C", "16C", "Blast")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(hjust = 0, vjust = 1, angle = 90),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.title = element_text(angle = 90),
        legend.position = c(-0.06, 1.13), 
        plot.margin = unit(c(3, 0.5, 0.5, 0.5), "cm")) +
  ggsave("bosTau7_LTR_log10rpm_expression.pdf") 

# plot as heatmap with fixed coordinate ratio
ggplot(data = rpm_df, aes(x = sample, y = repName)) + 
  geom_tile(aes(fill = rpm), colour = "white") + 
  scale_fill_gradient(low = "white", high = "black") +
  guides(fill = guide_colorbar(ticks = FALSE, 
                               title = expression(paste(log[10], " RPM")), 
                               title.position = "left", 
                               label.position = "left")) + 
  scale_x_discrete(position = "top", 
                   labels = c("GV", "MII", "4C", "8C", "16C", "Blast")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(hjust = 0, vjust = 1, angle = 90),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.title = element_text(angle = 90),
        legend.position = c(-0.4, 1.13), 
        plot.margin = unit(c(3, 0.5, 0.5, 0.5), "cm")) +
  coord_fixed(ratio = 1) + 
  ggsave("bosTau7_LTR_log10rpm_expression_fixed_ratio.pdf") 

  
