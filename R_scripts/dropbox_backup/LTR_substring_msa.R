# find signal (string) in sequences, substring them, realign and plot

# get name of sequence set
LTR_name <- "MTA_Mm"
LTR_set_name <- "expressed"
outfile <- paste0("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/", LTR_name, "_", LTR_set_name, ".fasta")

# get order in GV (3' half of LTR, FPKM)
FPKM_order <- fpkm_LTR$fullName[fpkm_LTR$fullName %in% names(LTR_seq_string)]

# get sequences, remove those without signal
LTR_seq <- as.character(readDNAStringSet(outfile, format = 'fasta'))
LTR_seq <- LTR_seq[str_detect(LTR_seq, "AATAAA")]

# locate string 
seq_string_location <- str_locate(LTR_seq, "AATAAA")
seq_string_location <- as.data.frame(seq_string_location)
seq_string_location$start <- seq_string_location$start - 20

# substring
LTR_seq_string <- sapply(X = 1:length(LTR_seq), FUN = function(X){
  substring(text = LTR_seq[X], first = seq_string_location[X, 1], last = seq_string_location[X, 2])
})

# write substring sequences as FASTA for align
write.fasta(as.list(LTR_seq_string),
            nbchar = 80,
            names = names(LTR_seq_string),
            as.string = TRUE,
            file.out = paste0(LTR_name, "_", LTR_set_name, "_polyA.fasta"),
            open = "w")

# aligne sequences by Clustal Omega
LTR_seq_string_aligned <- ClustalO(infile = paste0("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/", LTR_name, "_", LTR_set_name, "_polyA.fasta"), 
                                   outfile = paste0("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/", LTR_name, "_", LTR_set_name, "_polyA_aligned.fasta"), 
                                   threads = 8, 
                                   force = T)

# melt data.frame for plot
LTR_seq_string_aligned_df <- 
  as.data.frame(as.matrix(LTR_seq_string_aligned), stringsAsFactors = F) %>% 
  mutate(ID = rownames(.)) %>% 
  mutate(ID = factor(ID, rev(c("consensus", FPKM_order)))) %>% 
  dplyr::select(extract(., 1, ) %>% as.character() %>% equals("-") %>% not() %>% which()) %>%
  melt(id.vars = "ID") %>% 
  mutate(value = replace(value, value == "-", NA), 
         value = factor(value, levels = c("A", "C", "G", "T")))

# plot
ggplot(LTR_seq_string_aligned_df, aes(x = variable, y = ID)) +
  geom_tile(aes(fill = value)) +
  scale_fill_manual(values = color_pallete, 
                    breaks = c("A", "C", "G", "T")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 3, hjust = 1, vjust = 0),
        axis.ticks.y = element_blank()) +
  ggsave(paste0(LTR_name, "_", LTR_set_name, "_polyA_20nt_upstream_ordered.pdf"))
