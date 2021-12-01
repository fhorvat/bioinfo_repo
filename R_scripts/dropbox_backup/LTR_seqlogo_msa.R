# Hadleyverse
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

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/seq_logo/subclases")

################################################################################## functions
# Clustal Omega (vfranke)
ClustalO <- function(infile, outfile, ClustalO = NULL, threads = 12, 
                     what = 'DNA', force = T, format = 'fa'){
  
  if(!what %in% c('AA', 'DNA'))
    stop('can only align DNA or AA')
  if(is.null(ClustalO))
    ClustalO <- '/common/WORK/fhorvat/programi/clustal-omega-1.2.3/bin/clustalo'
  
  ### checks whether the file exists and whether to force the outfile
  ### if the file does exist and the force is off he reads the file
  if((file.exists(outfile) & force == T) | !file.exists(outfile)){
    cat('Running the alignmnent...\n')
    threads <- paste('--threads=', threads, sep = '')
    format <- paste('--outfmt=', format, sep = '')
    command <- paste(ClustalO, '-i', infile, '-o', outfile,'--force' ,threads, format)
    system(command)
  }
  
  cat('Returning the results...\n')
  if(what == 'AA')
    a <- readAAStringSet(outfile, format = 'fasta')
  if(what == 'DNA')
    a <- readDNAStringSet(outfile, format = 'fasta')
  
  return(a)
}

alignLTRClustalOAndPlot <- function(LTR_rptmsk, LTR_name){
  
  # get sequences 
  LTR_seq <- as.character(getSeq(x = Mmusculus, LTR_rptmsk))
  
  # write sequences as FASTA
  write.fasta(as.list(LTR_seq),
              nbchar = 80,
              names = LTR_rptmsk$fullName,
              as.string = TRUE,
              file.out = paste0(LTR_name, "_random_", length(LTR_rptmsk), ".fasta"),
              open = "w")
  
  # get sequences aligned by Clustal Omega
  LTR_seq_aligned <- ClustalO(infile = paste0(getwd(), "/", LTR_name, "_random_", length(LTR_rptmsk), ".fasta"), 
                              outfile = paste0(getwd(), "/", LTR_name, "_random_", length(LTR_rptmsk), "_aligned.fasta"), 
                              threads = 8, 
                              force = T)
  
  # melt data.frame for plot
  LTR_seq_aligned_df <- 
    as.data.frame(as.matrix(LTR_seq_aligned), stringsAsFactors = F) %>% 
    mutate(ID = rownames(.)) %>% 
    melt(id.vars = "ID") %>% 
    mutate(value = replace(value, value == "-", NA), 
           value = factor(value, levels = c("A", "C", "G", "T")))
  
  # plot
  ggplot(LTR_seq_aligned_df, aes(x = variable, y = ID)) +
    geom_tile(aes(fill = value)) +
    scale_fill_manual(values = color_pallete, 
                      breaks = c("A", "C", "G", "T")) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.title.y = element_blank(),
#           axis.text.y = element_text(size = 1, hjust = 1, vjust = 0),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggsave(paste0(LTR_name, "_random_", length(LTR_rptmsk), "_msa.pdf"))

  cat(paste0(LTR_name, " aligned and ploted\n\n"))
  
}

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

################################################################################## reading data
# mannualy set colors to the stages
# (A) green; (C) blue; (G) yellow; (T) red
color_pallete <- c("#66CD00", "#104E8B", "#FFD700", "#CD2626")
names(color_pallete) <- c("A", "C", "G", "T")

# LTR data
LTR_data <- 
  read_csv("LTR_ERVL_data_table.csv") %>% 
  dplyr::select(LTR_class, element_name = RepeatMasker_ID, consensus_width = consensus_sequence_length, consensus_length_fragment_numbers)

# repeatMasker
rptmsk <- read_delim("/common/WORK/fhorvat/reference/mm10/UCSC_repeatMasker_mm10_20161012.txt.gz", delim =  "\t")

# repeatMasker filtered
rptmsk_filter <- 
  right_join(rptmsk, LTR_data, by = "element_name") %>% 
  mutate(element_width = (end - start + 1)) %>% 
  filter(element_width > (consensus_width - (0.05 * consensus_width)), 
         element_width < (consensus_width + (0.05 * consensus_width))) %>% 
  dplyr::select(seqnames, start, end, strand, element_name, LTR_class) %>% 
  mutate(fullName = paste0(seqnames, ":",
                           start, "-", 
                           end, ":", 
                           strand, "|", 
                           element_name, "|", 
                           LTR_class))

# get fasta, align with ClustalO and plot (for each class)
for(LTR_class in unique(LTR_data$LTR_class)){
  
  set.seed(1234)
  LTR_rptmsk_sample <- 
    rptmsk_filter %>% 
    filter(LTR_class == LTR_class) %>% 
    sample_n(200) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
  
  alignLTRClustalOAndPlot(LTR_rptmsk = LTR_rptmsk_sample, LTR_name = LTR_class)
  
}

# plot polyA sequence logo
for(LTR_name in unique(LTR_data$LTR_class)){
  
  # read random 200 LTRs .fasta
  LTR_seq <- readDNAStringSet(filepath = paste0(getwd(), "/", LTR_name, "_random_200.fasta"), format = 'fasta')
  
  # get the most common signal start 
  signal_start <- 
    as.data.frame(vmatchPattern("AATAAA", LTR_seq)) %$%
    start %>% 
    remove_outliers() %>% 
    na.omit() %>% 
    median()
  
  # get sequence end
  seq_end <- ifelse(width(LTR_seq) < signal_start + 20, width(LTR_seq), signal_start + 20)
  
  # get the subsequence 40 nt around the signal start
  LTR_seq <- subseq(x = LTR_seq, start = signal_start - 20, end = seq_end)
  
  # write  .fasta 
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
  
  # get the most common signal start in aligned data 
  signal_start <- 
    as.data.frame(vmatchPattern("AATAAA", LTR_seq_string_aligned)) %$%
    start %>% 
    Mode() 
  
  # write logo with fixed signal start 
  weblogo(seqs = as.character(LTR_seq_string_aligned),
          open = F, 
          file.out = paste0(getwd(), "/", LTR_name, "_random_200", "_aligned_signal_logo.pdf"),
          format = "pdf", 
          color.scheme = "classic",
          lower = signal_start - 20, 
          upper = signal_start + 5, 
          show.xaxis = F, 
          errorbars = T)
  
  print(paste0(LTR_name, " logo done"))

}
