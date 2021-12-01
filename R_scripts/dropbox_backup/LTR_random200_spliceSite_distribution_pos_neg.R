################################################################################## reading data
# # repeatMasker
# rptmsk <-
#   read_delim("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_repeatMasker_mm10_20161012.txt.gz", delim =  "\t") %>%
#   dplyr::rename(LTR_subclass = element_name)
# 
# # LTR order
# LTR_order <- c("MLT1", "MLT2",
#                "ORR1G", "ORR1F", "ORR1E", "ORR1D", "ORR1C", "ORR1B", "ORR1A",
#                "MTE", "MTE2", "MTD", "MTC", "MTB", "MTA", "MT2A", "MT2B", "MT2C", "MT2")
# 
# # LTR data
# LTR_data <-
#   read_csv("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/LTR_alignment/seq_logo/LTR_ERVL_data_table.csv") %>%
#   dplyr::select(LTR_class, LTR_subclass = RepeatMasker_ID, consensus_width = consensus_sequence_length) %>%
#   mutate(LTR_class = replace(LTR_class, LTR_class == "MTB2", "MTB")) %>%
#   mutate(LTR_class = factor(LTR_class, levels = LTR_order))
# 
# # repeatMasker filtered
# rptmsk_filter <-
#   right_join(rptmsk, LTR_data, by = "LTR_subclass") %>%
#   mutate(element_width = (end - start + 1)) %>%
#   filter(element_width > (consensus_width - (0.05 * consensus_width)),
#          element_width < (consensus_width + (0.05 * consensus_width))) %>%
#   dplyr::select(seqnames, start, end, strand, LTR_subclass, LTR_class) %>%
#   mutate(fullName = paste0(seqnames, ":",
#                            start, "-",
#                            end, ":",
#                            strand, "|",
#                            LTR_subclass, "|",
#                            LTR_class))
# 
# #get all LTRs from one class
# LTR_random <-
#   makeGRangesFromDataFrame(rptmsk_filter, keep.extra.columns = T)

################################################################################## 
# get 200 random LTRs from one class
LTR_random <-
  read_csv("LTR_random200_perClass.csv") %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

# get sequences 
LTR_seqences <- getSeq(x = Mmusculus, LTR_random)
names(LTR_seqences) <- LTR_random$fullName

############################################################### negative
distribution_all <-
  unlist(vmatchPattern(splice_pattern, LTR_seqences, fixed = F)) %>% 
  as.data.frame() %>% 
  dplyr::select(start, names) %>% 
  separate(names, c("LTR_range", "LTR_subclass", "LTR_class"), sep = "\\|")  %>% 
  separate(LTR_range, c("LTR_seqnames", "LTR_range", "LTR_strand"), sep = ":") %>% 
  separate(LTR_range, c("LTR_start", "LTR_end"), sep = "-") %>% 
  mutate(LTR_start = as.numeric(LTR_start), 
         LTR_end = as.numeric(LTR_end), 
         LTR_width = (LTR_end - LTR_start) + 1, 
         pattern = splice_pattern) %>% 
  left_join(LTR_data, by = "LTR_subclass") %>% 
  dplyr::select(start, LTR_subclass, LTR_class = LTR_class.x, LTR_width, consensus_width, LTR_strand, pattern) %>% 
  mutate(start = -(LTR_width - start),
         LTR_class = factor(LTR_class, levels = LTR_order)) %>% 
  dplyr::filter(start >= -(consensus_width)) 

consensus_length_negative <- 
  LTR_data %>% 
  group_by(LTR_class) %>% 
  summarise(consensus_width = max(consensus_width)) %>% 
  mutate(consensus_width = -consensus_width)

# get splice start median per class
signal_start <- 
  distribution_all %>%
  group_by(LTR_class) %>% 
  summarize(start = round(median(start)))

# plot pattern start distribution as jitterplot
ggplot(data = distribution_all, aes(x = LTR_class, start)) +
  geom_bar(data = consensus_length_negative, aes(x = LTR_class, y = consensus_width), fill = "#D9D9D9", stat = "identity") +
  geom_boxplot(fill = NA, show.legend = F, outlier.shape = NA) +
  geom_jitter(aes(colour = pattern),  size = 0.8, width = 0.25) +
  geom_text(data = signal_start, aes(x = LTR_class, y = start, label = start), colour = "white", size = 3, fontface = "bold") +
  xlab(NULL) +
  ylab("length(nt)") +
  ylim(c(min(consensus_length_negative$consensus_width), 0)) +
  scale_colour_manual(values = "red") +
  scale_x_discrete(drop = FALSE, position = "bottom") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 2))) + 
  ggsave(filename = paste0("LTR_random200_", splice_pattern, "_plot_negative.png"), height = 15, width = 30, unit = "cm", dpi = 500)


############################################################### positive
distribution_all <-
  unlist(vmatchPattern(splice_pattern, LTR_seqences, fixed = F)) %>% 
  as.data.frame() %>% 
  dplyr::select(start, names) %>% 
  separate(names, c("LTR_range", "LTR_subclass", "LTR_class"), sep = "\\|")  %>% 
  separate(LTR_range, c("LTR_seqnames", "LTR_range", "LTR_strand"), sep = ":") %>% 
  separate(LTR_range, c("LTR_start", "LTR_end"), sep = "-") %>% 
  mutate(LTR_start = as.numeric(LTR_start), 
         LTR_end = as.numeric(LTR_end), 
         LTR_width = (LTR_end - LTR_start) + 1, 
         pattern = splice_pattern) %>% 
  left_join(LTR_data, by = "LTR_subclass") %>% 
  dplyr::select(start, LTR_subclass, LTR_class = LTR_class.x, LTR_width, consensus_width, LTR_strand, pattern) %>% 
  mutate(start = -(LTR_width - start),
         LTR_class = factor(LTR_class, levels = LTR_order)) %>% 
  dplyr::filter(start >= -(consensus_width)) %>% 
  mutate(start = LTR_width + start)

consensus_length_negative <- 
  LTR_data %>% 
  group_by(LTR_class) %>% 
  summarise(consensus_width = max(consensus_width)) 

# get splice start median per class
signal_start <- 
  distribution_all %>%
  group_by(LTR_class) %>% 
  summarize(start = round(median(start)))

# plot pattern start distribution as jitterplot
ggplot(data = distribution_all, aes(x = LTR_class, start)) +
  geom_bar(data = consensus_length_negative, aes(x = LTR_class, y = consensus_width), fill = "#D9D9D9", stat = "identity") +
  geom_boxplot(fill = NA, show.legend = F, outlier.shape = NA) +
  geom_jitter(aes(colour = pattern),  size = 0.8, width = 0.25) +
  geom_text(data = signal_start, aes(x = LTR_class, y = start, label = start), colour = "white", size = 3, fontface = "bold") +
  xlab(NULL) +
  ylab("length(nt)") +
  ylim(c(0, max(consensus_length_negative$consensus_width))) +
  scale_colour_manual(values = "red") +
  scale_x_discrete(drop = FALSE, position = "bottom") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 2))) + 
  ggsave(filename = paste0("LTR_random200_", splice_pattern, "_plot_positive.png"), height = 15, width = 30, unit = "cm", dpi = 500)
