# RWebLogo
library(RWebLogo)
library(seqinr)

LTR_seq <- readDNAStringSet(filepath = paste0(getwd(), "/", LTR_name, "_random_200.fasta"), format = 'fasta')

signal_start <- 
  as.data.frame(vmatchPattern("AATAAA", LTR_seq)) %$%
  start %>% 
  table() %>% 
  sort(decreasing = T) %>% 
  .[1] %>% 
  names() %>% 
  as.integer()

LTR_seq <- subseq(x = LTR_seq, start = signal_start - 20, end = signal_start + 10)

write.fasta(as.list(LTR_seq),
            nbchar = 80,
            names = names(LTR_seq),
            as.string = TRUE,
            file.out = paste0(LTR_name, "_random_200_polyA.fasta"),
            open = "w")

# align sequences by Clustal Omega
LTR_seq_string_aligned <- ClustalO(infile = paste0(getwd(), "/", LTR_name, "_random_200_polyA.fasta"), 
                                   outfile = paste0(getwd(), "/", LTR_name, "_random_200_polyA_aligned.fasta"), 
                                   threads = 8, 
                                   force = T)

signal_start <- 
  as.data.frame(vmatchPattern("AATAAA", LTR_seq_string_aligned)) %$%
  start %>% 
  table() %>% 
  sort(decreasing = T) %>% 
  .[1] %>% 
  names() %>% 
  as.integer()

weblogo(seqs = as.character(LTR_seq_string_aligned),
        open = F, 
        file.out = paste0(getwd(), "/", LTR_name, "_random_200", "_aligned_signal_logo.png"),
        format = "png", 
        color.scheme = "classic", 
        lower = signal_start - 20, 
        upper = signal_start + 5, 
        show.xaxis = F)

