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
changeWD("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/LTR_alignment/polyA_splice_pattern_distribution/random_LTRs/test")

# finds mode of the vector (the most common value)
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Clustal Omega (vfranke)
ClustalO <- function(infile, outfile, ClustalO = NULL, threads = 12, 
                     what = 'DNA', force = T, format = 'fa', iterations = 0){
  
  if(!what %in% c('AA', 'DNA'))
    stop('can only align DNA or AA')
  if(is.null(ClustalO))
    ClustalO <- '/common/WORK/fhorvat/programi/clustal-omega-1.2.3/bin/clustalo'
  
  ### checks whether the file exists and whether to force the outfile
  ### if the file does exist and the force is off he reads the file
  if((file.exists(outfile) & force == T) | !file.exists(outfile)){
    cat(gsub(".*\\/", "", infile), "running the alignmnent", "\n")
    threads <- paste('--threads=', threads, sep = '')
    format <- paste('--outfmt=', format, sep = '')
    iterations <- paste("--iter=", iterations, sep ="")
    pileup <- paste("--pileup", sep = "")
    command <- paste(ClustalO, '-i', infile, '-o', outfile,'--force', threads, format, iterations)
    system(command)
  }
  
  cat(gsub(".*\\/", "", infile), "alignment done", "\n")
  if(what == 'AA')
    a <- readAAStringSet(outfile, format = 'fasta')
  if(what == 'DNA')
    a <- readDNAStringSet(outfile, format = 'fasta')
  
  return(a)
}

# plots aligned sequences as "heatmap"
msaPlot <- function(aligned_sequences, plot_file_name){
  
  # melt data.frame for plot
  aligned_sequences_df <- 
    as.data.frame(as.matrix(aligned_sequences), stringsAsFactors = F) %>% 
    mutate(ID = rownames(.)) %>% 
    melt(id.vars = "ID") %>% 
    mutate(value = replace(value, value == "-", NA), 
           value = factor(value, levels = c("A", "C", "G", "T")))
  
  # plot
  ggplot(aligned_sequences_df, aes(x = variable, y = ID)) +
    geom_tile(aes(fill = value)) +
    scale_fill_manual(values = color_pallete, 
                      breaks = c("A", "C", "G", "T")) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 3, hjust = 1, vjust = 0, angle = 90),
          axis.ticks.x = element_blank(), 
          axis.title.y = element_blank(),
          #           axis.text.y = element_text(size = 1, hjust = 1, vjust = 0),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggsave(plot_file_name)
  
  cat(paste0(plot_file_name, " done\n"))
  
}

# mannualy set colors to the stages
# (A) green; (C) blue; (G) yellow; (T) red
color_pallete <- 
  c("#66CD00", "#104E8B", "#FFD700", "#CD2626") %>% 
  set_names(c("A", "C", "G", "T"))

################################################################################## reading data
# LTR order 
LTR_order <- c("MLT1", "MLT2", 
               "ORR1G", "ORR1F", "ORR1E", "ORR1D", "ORR1C", "ORR1B", "ORR1A", 
               "MTE", "MTE2", "MTD", "MTC", "MTB", "MTA", "MT2A", "MT2B", "MT2C", "MT2")

# LTR data
LTR_data <- 
  read_csv("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/LTR_alignment/seq_logo/LTR_ERVL_data_table.csv") %>% 
  dplyr::select(LTR_class, LTR_subclass = RepeatMasker_ID, consensus_width = consensus_sequence_length) %>% 
  mutate(LTR_class = replace(LTR_class, LTR_class == "MTB2", "MTB")) %>% 
  mutate(LTR_class = factor(LTR_class, levels = LTR_order)) 

################################################################################## 
# get 200 random LTRs from one class
LTR_random <- 
  read_csv("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/LTR_alignment/data_tables/LTR_random200_perClass.csv") %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

# get sequences 
LTR_seqences <- 
  getSeq(x = Mmusculus, LTR_random) %>% 
  set_names(LTR_random$fullName)

################################################################################## get pattern start positions
# splice pattern
splice_pattern <- "TGTAAGY"

# get splice pattern start positions - 3' end
splice_distribution <- 
  unlist(vmatchPattern(splice_pattern, LTR_seqences, fixed = F)) %>% 
  as.data.frame() %>% 
  dplyr::select(start, names) %>% 
  separate(names, c("LTR_range", "LTR_subclass", "LTR_class"), sep = "\\|") %>% 
  separate(LTR_range, c("LTR_seqnames", "LTR_range", "LTR_strand"), sep = ":") %>% 
  separate(LTR_range, c("LTR_start", "LTR_end"), sep = "-") %>% 
  mutate(LTR_start = as.numeric(LTR_start), 
         LTR_end = as.numeric(LTR_end), 
         LTR_width = (LTR_end - LTR_start) + 1) %>% 
  left_join(LTR_data, by = "LTR_subclass") %>% 
  dplyr::select(start, LTR_subclass, LTR_class = LTR_class.x, LTR_width, consensus_width, LTR_strand) %>% 
  mutate(start = -(LTR_width - start),
         LTR_class = factor(LTR_class, levels = LTR_order)) %>%
  dplyr::filter(start >= -(consensus_width)) %>% 
  mutate(start = start + LTR_width)
# %>% dplyr::count(LTR_class)

# get splice start median per class
signal_start <- 
  splice_distribution %>%
  group_by(LTR_class) %>% 
  summarize(start = round(median(start)))

# sequence width
seq_width <- 10

# names of LTRs with splice pattern 
splice_distribution_names <- 
  unlist(vmatchPattern(splice_pattern, LTR_seqences, fixed = F)) %>% 
  names()

# loop over LTR classes
for(file_name in LTR_order){
  
  # get signal start
  signal_start_per_class <- 
    signal_start %>% 
    dplyr::filter(LTR_class == file_name) %$%
    start
  
  # subset sequences by LTR class
  LTR_sequences_per_class <- LTR_seqences[grep(paste0("\\|", file_name, "$"), names(LTR_seqences))]
  
  # order so that the ones with splice site pattern go first
  LTR_sequences_per_class <- c(LTR_sequences_per_class[names(LTR_sequences_per_class) %in% splice_distribution_names],
                               LTR_sequences_per_class[!names(LTR_sequences_per_class) %in% splice_distribution_names])

  # get sequences end (either signal start + 20 or sequence end)
  sub_sequences_end <- ifelse(width(LTR_sequences_per_class) < signal_start_per_class + seq_width, width(LTR_sequences_per_class), signal_start_per_class + seq_width)
  
  # get the subsequence 20 nt around the signal start
  sub_sequences <- subseq(x = LTR_sequences_per_class, start = signal_start_per_class - seq_width, end = sub_sequences_end)
  
  # write .fasta
  write.fasta(as.list(sub_sequences),
              nbchar = 80,
              names = names(sub_sequences),
              as.string = TRUE,
              file.out = paste0(file_name, "_", seq_width, ".fasta"),
              open = "w")

  # align sequences by Clustal Omega
  sub_sequences_aligned <- ClustalO(infile = paste0(getwd(), "/", file_name, "_", seq_width, ".fasta"),
                                    outfile = paste0(getwd(), "/", file_name, "_", seq_width,"_msa_iter.fasta"),
                                    threads = 8,
                                    force = T,
                                    iterations = 100)
  
  # plot alignment
  msaPlot(aligned_sequences = sub_sequences_aligned, 
          plot_file_name = paste0(file_name, "_", seq_width, "_msa_iter_plot.pdf"))
  
  # get gap index
  gap_index <-
    consensusMatrix(sub_sequences_aligned, as.prob = F, baseOnly = T)[5, ] %>%
    is_greater_than(190) %>%
    which()

  if(length(gap_index) > 0){

    # keep only those columns in alignment
    sub_sequences_aligned <-
      as.matrix(sub_sequences_aligned) %>%
      .[, -gap_index] %>%
      apply(., 1, paste, collapse = "") %>%
      DNAStringSet()

  }
  
  # find the most common pattern start in aligned sequences
  signal_start_aligned <- 
    as.data.frame(vmatchPattern(splice_pattern, sub_sequences_aligned, max.mismatch = 1)) %$%
    start %>% 
    Mode() 
  
  if(!is.na(signal_start_aligned)){
    
    if(signal_start_aligned < 11){
      signal_start_aligned <- 11
    }
    
    if(signal_start_aligned > (max(width(sub_sequences_aligned)) - 10)){
      signal_start_aligned <- max(width(sub_sequences_aligned)) - 10
    }
    
    sub_sequences_aligned <- subseq(x = sub_sequences_aligned, start = signal_start_aligned - seq_width, end = signal_start_aligned + seq_width)
    
  }
  
  # plot alignment short
  msaPlot(aligned_sequences = sub_sequences_aligned,
          plot_file_name = paste0(file_name, "_", seq_width, "_msa_iter_plot_short.pdf"))
  
  # write logo with fixed signal start 
  weblogo(seqs = as.character(sub_sequences_aligned),
          open = F, 
          file.out = paste0(getwd(), "/", file_name, "_", seq_width, "_msa_iter_logo.pdf"),
          format = "pdf", 
          color.scheme = "classic",
          show.xaxis = F, 
          errorbars = F)
}



