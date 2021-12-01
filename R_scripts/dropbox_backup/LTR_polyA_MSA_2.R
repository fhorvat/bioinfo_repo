# get aligned file
LTR_name <- "MTA_Mm"
LTR_set_name <- "expressed"
outfile <- paste0("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/", LTR_name, "_", LTR_set_name, "_aligned.fasta")

# get sequences 
LTR_seq <- readDNAStringSet(outfile, format = 'fasta')

LTR_seq_char <- as.character(LTR_seq)
LTR_seq_char <- LTR_seq_char[str_detect(LTR_seq_char, "AATAAA")]

LTR_seq <- LTR_seq[names(LTR_seq) %in% names(LTR_seq_char)]
FPKM_order <- fpkm_LTR$fullName[fpkm_LTR$fullName %in% names(LTR_seq)]

# locate string
polyA_signal <- str_locate(LTR_seq_char[1], "AATAAA")
polyA_signal <- as.data.frame(polyA_signal)

polyA_signal_shift <- polyA_signal
polyA_signal_shift$start <- polyA_signal_shift$start - 20

# convert to DNAStringSet and write as fasta
LTR_seq_polyA <- 
  as.data.frame(as.matrix(LTR_seq), stringsAsFactors = F) %>% 
  dplyr::select(polyA_signal_shift$start:polyA_signal_shift$end) 
LTR_seq_polyA <- do.call("paste0", LTR_seq_polyA)
names(LTR_seq_polyA) <- names(LTR_seq)
LTR_seq_polyA <- LTR_seq_polyA[match(FPKM_order, names(LTR_seq_polyA))]

# # write as fasta
# LTR_seq_polyA <- DNAStringSet(LTR_seq_polyA)
# write.fasta(as.list(as.character(LTR_seq_polyA)),
#             nbchar = 80,
#             names = names(LTR_seq_polyA),
#             as.string = TRUE,
#             file.out = paste0(LTR_name, "_", LTR_set_name, "_polyA_upstream_with_consensus2.fasta"),
#             open = "w")

# write as text
writeLines(paste0(LTR_seq_polyA, "    ", names(LTR_seq_polyA)), 
           paste0(LTR_name, "_", LTR_set_name, "_polyA_upstream_with_consensus_ordered.txt"))

# melt data.frame for plot
LTR_seq_polyA_df <- 
  as.data.frame(as.matrix(LTR_seq), stringsAsFactors = F) %>% 
  dplyr::select(polyA_signal_shift$start:polyA_signal_shift$end) %>% 
  set_colnames(polyA_signal_shift$start:polyA_signal_shift$end) %>% 
  mutate(ID = rownames(.)) %>% 
  mutate(ID = factor(ID, rev(c("consensus", FPKM_order)))) %>% 
  dplyr::select(extract(., 1, ) %>% as.character() %>% equals("-") %>% not() %>% which()) %>% 
  melt(id.vars = "ID") %>% 
  mutate(value = replace(value, value == "-", NA), 
         value = factor(value, levels = c("A", "C", "G", "T")))

# plot
ggplot(LTR_seq_polyA_df, aes(x = variable, y = ID)) +
  geom_tile(aes(fill = value)) +
  scale_fill_manual(values = color_pallete, 
                    breaks = c("A", "C", "G", "T")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 3, hjust = 1, vjust = 0),
        axis.ticks.y = element_blank()) +
  ggsave(paste0(LTR_name, "_", LTR_set_name, "_polyA_upstream_with_consensus_ordered.pdf"))
