LTR_name <- "MTA_Mm"
LTR_set_name <- "expressed"
outfile <- paste0("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/", LTR_name, "_", LTR_set_name, "_aligned.fasta")

# get sequences 
LTR_seq_original <- as.character(readDNAStringSet(outfile, format = 'fasta'))
LTR_seq <- str_replace_all(LTR_seq_original, "-", "")
names(LTR_seq) <- names(LTR_seq_original)
LTR_seq <- LTR_seq[str_detect(LTR_seq, "AATAAA")]

# locate string
polyA_signal <- str_locate(LTR_seq, "AATAAA")
polyA_signal <- as.data.frame(polyA_signal)
polyA_signal_shift <- polyA_signal
polyA_signal_shift$start <- polyA_signal_shift$start - 20

# substring
LTR_seq_polyA <- sapply(X = 1:length(LTR_seq), FUN = function(X){
  substring(text = LTR_seq[X], first = polyA_signal_shift[X, 1], last = polyA_signal_shift[X, 2])
})

LTR_seq_polyA <- DNAStringSet(LTR_seq_polyA)

# get sequences 
LTR_seq_polyA <- as.character(LTR_seq_polyA)

# write sequences as FASTA
write.fasta(as.list(LTR_seq_polyA),
            nbchar = 80,
            names = names(LTR_seq_polyA),
            as.string = TRUE,
            file.out = paste0(LTR_name, "_", LTR_set_name, "_polyA.fasta"),
            open = "w")

# get sequences aligned by Clustal Omega
LTR_seq_polyA_aligned <- ClustalO(infile = paste0("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/", LTR_name, "_", LTR_set_name, "_polyA.fasta"), 
                                  outfile = paste0("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/", LTR_name, "_", LTR_set_name, "_polyA_aligned.fasta"), 
                                  threads = 8, 
                                  force = T)

# order by FPKM in GV
FPKM_order <- fpkm_LTR$fullName[fpkm_LTR$fullName %in% names(LTR_seq_polyA)]

# melt data.frame for plot
LTR_seq_polyA_df <- 
  as.data.frame(as.matrix(LTR_seq_polyA_aligned), stringsAsFactors = F) %>% 
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
  ggsave(paste0(LTR_name, "_", LTR_set_name, "_polyA_upstream_with_consensus_ordered_4.pdf"))
