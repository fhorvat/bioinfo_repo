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

changeWD("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/seq_logo/pattern_distribution")

################################################################################## reading data
# LTR data
LTR_data <- 
  read_csv("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/seq_logo/LTR_ERVL_data_table.csv") %>% 
  dplyr::select(LTR_class, element_name = RepeatMasker_ID, consensus_width = consensus_sequence_length, consensus_length_fragment_numbers)

# repeatMasker
rptmsk <- read_delim("/common/WORK/fhorvat/reference/mm10/UCSC_repeatMasker_mm10_20161012.txt.gz", delim =  "\t")

# repeatMasker filtered
rptmsk_filter <- 
  right_join(rptmsk, LTR_data, by = "element_name") %>% 
  mutate(element_width = (end - start + 1)) %>% 
  filter(element_width > (consensus_width - (0.05 * consensus_width)), 
         element_width < (consensus_width + (0.05 * consensus_width))) %>% 
  dplyr::select(seqnames, start, end, strand, element_name, LTR_class, consensus_width) %>% 
  mutate(fullName = paste0(seqnames, ":",
                           start, "-", 
                           end, ":", 
                           strand, "|", 
                           element_name, "|", 
                           LTR_class))

# get consensus length for each class
consenus_length <-
  rptmsk_filter %>%
  group_by(LTR_class) %>%
  summarize(consensus_width = max(consensus_width)) %>% 
  mutate(consensus_width = ceiling(consensus_width + (0.05 * consensus_width)))

################################################################################## reading data
# get x random LTRs from one class
sample_size <- "all"

set.seed(1234)
LTR_rptmsk_sample <- 
  rptmsk_filter %>% 
  group_by(LTR_class) %>% 
  mutate(group_size = n()) %>% 
  # filter(fullName %in% sample(fullName, ifelse(group_size > sample_size, sample_size, group_size))) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

# get sequences 
LTR_seqences <- 
  getSeq(x = Mmusculus, LTR_rptmsk_sample) %>% 
  setNames(LTR_rptmsk_sample$fullName)

################################################################################## set pattern
# strict splice site 
pattern <- "AAGTGTAA"

# relaxed splice site (A-G/C-T-G-T-A/G-A) written as IUPAC ambiguous DNA
pattern <- "ASTGTRA"

################################################################################## get pattern distribution
# get sequences with signal
signal_start <- 
  unlist(vmatchPattern(pattern, LTR_seqences, fixed = F)) %>% 
  as.data.frame() %>%
  separate(names, c("LTR_range", "LTR_subclass", "LTR_class"), sep = "\\|") %>% 
  separate(LTR_range, c("seqnames", "range", "strand"), sep = ":") %>% 
  separate(range, c("LTR_start", "LTR_end"), sep = "-") %>% 
  mutate(LTR_class = factor(LTR_class, levels = unique(rptmsk_filter$LTR_class)))

################################################################################## plot
# plot pattern start distribution as jitterplot
ggplot(data = signal_start, aes(LTR_class, start)) +
  geom_bar(data = consenus_length, aes(x = factor(LTR_class), y = consensus_width), stat = "identity", alpha = 0.2) +
  geom_jitter(aes(colour = LTR_class), size = 0.1, width = 0.25, show.legend = F) +
  # stat_summary(fun.y = "median", colour = "red", size = 1, geom = "point") +
  # geom_text(data = signal_start_med, aes(x = factor(LTR_class), y = start_med, label = start_med), size = 3, vjust = 1.5) +
  xlab("LTR class") +
  ylab(paste0(pattern, " start position")) +
  ylim(c(0, max(consenus_length$consensus_width))) +
  scale_colour_discrete(drop = FALSE) +
  scale_x_discrete(drop = FALSE) + 
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0)) +
  ggsave(paste0("spliceSite_distribution_", sample_size, "_test.png"))

##################################################################################
##################################################################################
# # plot pattern start distribution as boxplot
# ggplot(data = signal_start, aes(factor(LTR_class), start)) +
#   geom_boxplot(aes(colour = LTR_class), show.legend = F) +
#   geom_text(data = signal_start_med, aes(x = factor(LTR_class), y = start_med, label = start_med), size = 3) + 
#   xlab("LTR class") +
#   ylab("AATAAA start position") +
#   theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0)) +
#   ggsave(paste0("class_signal_start_boxplot_absolute_", sample_size, ".pdf"))

# # calculate median for each class
# signal_start_med <- 
#   signal_start %>% 
#   group_by(LTR_class) %>% 
#   summarize(start_med = round(median(start)))

# # stringr version of pattern matching
# pattern <- "AA[G,C]TGT[A,G]A"
# signal_start <-
#   str_locate_all(LTR_seqences_chr, pattern) %>%
#   setNames(names(LTR_seqences_chr)) %>%
#   lapply(., as.data.frame, stringAsFactors = F) %>%
#   bind_rows(., .id = "names") %>%
#   dplyr::select(start, end, names) %>%
#   separate(names, c("LTR_range", "LTR_subclass", "LTR_class"), sep = "\\|") %>%
#   separate(LTR_range, c("seqnames", "range", "strand"), sep = ":") %>%
#   separate(range, c("LTR_start", "LTR_end"), sep = "-")
