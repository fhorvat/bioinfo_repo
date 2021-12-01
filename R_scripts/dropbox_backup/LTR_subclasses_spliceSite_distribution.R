# Hadley-verse
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(reshape2)
library(magrittr)

# Bioconducor libraries
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(BiocParallel)
library(BSgenome.Mmusculus.UCSC.mm10)
library(seqinr)
library(Biostrings)
library(RWebLogo)

options(bitmapType = "cairo")

changeWD <- function(wd){
  if(!dir.exists(wd)){
    dir.create(wd, showWarnings = TRUE, recursive = FALSE, mode = "0755")
    setwd(wd)
    cat("Working directory created and set to", getwd(), "\n")
  } else {
    setwd(wd)
    cat("Working directory set to", getwd(), "\n")
  }
}

changeWD("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/seq_logo/random200_LTR/")

################################################################################## reading data
# repeatMasker
rptmsk <- 
  read_delim("/common/WORK/fhorvat/reference/mm10/UCSC_repeatMasker_mm10_20161012.txt.gz", delim =  "\t") %>% 
  dplyr::rename(LTR_subclass = element_name)

# LTR order 
LTR_order <- c("MLT1", "MLT2", 
               "ORR1G", "ORR1F", "ORR1E", "ORR1D", "ORR1C", "ORR1B", "ORR1A", 
               "MTE", "MTE2", "MTD", "MTC", "MTB", "MTA", "MT2A", "MT2B", "MT2C", "MT2")

# LTR data 
LTR_data <- 
  read_csv("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/seq_logo/LTR_ERVL_data_table.csv") %>% 
  dplyr::select(LTR_class, LTR_subclass = RepeatMasker_ID, consensus_width = consensus_sequence_length, consensus_length_fragment_numbers) %>% 
  mutate(LTR_class = replace(LTR_class, LTR_class == "MTB2", "MTB")) 

# repeatMasker filtered
rptmsk_filter <- 
  right_join(rptmsk, LTR_data, by = "LTR_subclass") %>% 
  mutate(element_width = (end - start + 1)) %>% 
  filter(element_width > (consensus_width - (0.05 * consensus_width)), 
         element_width < (consensus_width + (0.05 * consensus_width))) %>% 
  dplyr::select(seqnames, start, end, strand, LTR_subclass, LTR_class) %>% 
  mutate(fullName = paste0(seqnames, ":",
                           start, "-", 
                           end, ":", 
                           strand, "|", 
                           LTR_subclass, "|", 
                           LTR_class))

################################################################################## 
# get 200 random LTRs from one class
LTR_random <- makeGRangesFromDataFrame(rptmsk_filter, keep.extra.columns = T)
# write_csv(x = as.data.frame(LTR_random), path = "LTR_random200_perClass.csv")

# get sequences 
LTR_seqences <- getSeq(x = Mmusculus, LTR_random)
names(LTR_seqences) <- LTR_random$fullName

################################################################################## get polyA pattern distribution
# set pattern
seq_pattern <- "AANTGTAAG"
# seq_pattern <- "KGTRAG"

# get splice pattern start positions (T/G-G-T-A/G-A-G) 
splice_distribution <- 
  unlist(vmatchPattern(seq_pattern, LTR_seqences, fixed = F)) %>% 
  as.data.frame() %>%
  dplyr::select(start, names) %>% 
  separate(names, c("LTR_range", "LTR_subclass", "LTR_class"), sep = "\\|") %>% 
  left_join(LTR_data, by = "LTR_subclass") %>% 
  dplyr::filter(start < consensus_width) %>% 
  mutate(pattern = seq_pattern)

# get polyA pattern start positions
polyA_distribution <- 
  unlist(vmatchPattern("AATAAA", LTR_seqences, fixed = F)) %>% 
  as.data.frame() %>%
  dplyr::select(start, names) %>% 
  separate(names, c("LTR_range", "LTR_subclass", "LTR_class"), sep = "\\|") %>% 
  left_join(LTR_data, by = "LTR_subclass") %>% 
  dplyr::filter(start < consensus_width) %>% 
  mutate(pattern = "polyA")

################################################################################## plot
# join polyA/splice distribution
distribution_all <-
  rbind(polyA_distribution, splice_distribution) %>% 
  mutate(start = -(consensus_width - start), 
         pattern = factor(pattern, levels = c(seq_pattern, "polyA")), 
         LTR_subclass = factor(LTR_subclass, levels = unique(LTR_subclass)))

# LTR data with consensus length (Maja)
consensus_length_negative <- 
  LTR_data %>% 
  mutate(consensus_width = -consensus_width) %>% 
  dplyr::select(LTR_subclass, consensus_width) %>% 
  mutate(LTR_subclass = factor(LTR_subclass, levels = levels(distribution_all$LTR_subclass)))

# plot pattern start distribution as jitterplot
ggplot(data = distribution_all, aes(x = LTR_subclass, start)) +
  geom_bar(data = consensus_length_negative, aes(x = LTR_subclass, y = consensus_width), stat = "identity", alpha = 0.2) +
  geom_jitter(aes(colour = pattern), size = 0.1, width = 0.25) +
  xlab("LTR class") +
  ylab("pattern start position") +
  ylim(c(min(consensus_length_negative$consensus_width), 0)) +
  scale_colour_manual(values = c("red", "black")) +
  scale_x_discrete(drop = T) + 
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggsave(paste0("LTR_subclasses_", seq_pattern, "_distribution_plot.png"), width = 14, height = 8)
