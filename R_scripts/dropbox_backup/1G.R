rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd(".")

######################################################## LIBRARIES
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(reshape2)
library(magrittr)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(BiocParallel)
library(BSgenome.Mmusculus.UCSC.mm10)
library(seqinr)
library(Biostrings)

######################################################## FUNCTIONS

######################################################## READ DATA
# LTR order 
LTR_order <- c("MLT1", "MLT2", 
               "ORR1G", "ORR1F", "ORR1E", "ORR1D", "ORR1C", "ORR1B", "ORR1A", 
               "MTE", "MTE2", "MTD", "MTC", "MTB", "MTA", "MT2A", "MT2B", "MT2C", "MT2")

# LTR data (all LTR classes, subclasses and consensus length for each element)
LTR_data <- 
  read_csv("LTR_ERVL_data_table.csv") %>% 
  dplyr::select(LTR_class, LTR_subclass = RepeatMasker_ID, consensus_width = consensus_sequence_length) %>% 
  mutate(LTR_class = replace(LTR_class, LTR_class == "MTB2", "MTB")) %>% 
  mutate(LTR_class = factor(LTR_class, levels = LTR_order)) 

# coordinates of 200 random LTRs for each class
LTR_random <-
  read_csv("LTR_random200_perClass.csv") %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

######################################################## MAIN CODE
# sequences of 200 random LTRs
LTR_seqences <- getSeq(x = Mmusculus, LTR_random)
names(LTR_seqences) <- LTR_random$fullName

# splice donor motif start positions - RED
splice_distribution <- 
  unlist(vmatchPattern("TGTAAGY", LTR_seqences, fixed = F)) %>% 
  as.data.frame() %>% 
  dplyr::select(start, names) %>% 
  separate(names, c("LTR_range", "LTR_subclass", "LTR_class"), sep = "\\|") %>% 
  mutate(pattern = "TGTAAGY")

# polyA motif start positions - GREY
polyA_distribution <- 
  unlist(vmatchPattern("AATAAA", LTR_seqences, fixed = F)) %>% 
  as.data.frame() %>%
  dplyr::select(start, names) %>% 
  separate(names, c("LTR_range", "LTR_subclass", "LTR_class"), sep = "\\|") %>% 
  mutate(pattern = "polyA")

# both motif start positions together, set 3' end as 0 position
distribution_all <-
  rbind(polyA_distribution, splice_distribution) %>% 
  separate(LTR_range, c("LTR_seqnames", "LTR_range", "LTR_strand"), sep = ":") %>% 
  separate(LTR_range, c("LTR_start", "LTR_end"), sep = "-") %>% 
  mutate(LTR_start = as.numeric(LTR_start), 
         LTR_end = as.numeric(LTR_end), 
         LTR_width = (LTR_end - LTR_start) + 1) %>% 
  left_join(LTR_data, by = "LTR_subclass") %>% 
  dplyr::select(start, LTR_subclass, LTR_class = LTR_class.x, LTR_width, consensus_width, LTR_strand, pattern) %>% 
  mutate(start = -(LTR_width - start), 
         pattern = factor(pattern, levels = c("polyA", "TGTAAGY")),
         LTR_class = factor(LTR_class, levels = LTR_order)) %>% 
  dplyr::filter(start >= -(consensus_width))

# consensus length for each class
consensus_length_negative <- 
  LTR_data %>% 
  group_by(LTR_class) %>% 
  summarise(consensus_width = max(consensus_width)) %>% 
  mutate(consensus_width = -consensus_width)

# plot pattern start distribution
ggplot(data = distribution_all, aes(x = LTR_class, start)) +
  geom_bar(data = consensus_length_negative, aes(x = LTR_class, y = consensus_width), fill = "#D9D9D9", stat = "identity") +
  geom_jitter(aes(colour = pattern),  size = 0.8, width = 0.25) +
  xlab(NULL) +
  ylab("length(nt)") +
  ylim(c(min(consensus_length_negative$consensus_width), 0)) +
  scale_colour_manual(values = c("grey40", "red")) +
  scale_x_discrete(drop = FALSE, position = "bottom") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 2))) + 
  ggsave(filename = "LTR_random200_splice_pattern_plot.png", height = 15, width = 30, unit = "cm", dpi = 500)
