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

changeWD("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/seq_logo/test")

################################################################################## functions
# mannualy set colors to the stages
# (A) green; (C) blue; (G) yellow; (T) red
color_pallete <- c("#66CD00", "#104E8B", "#FFD700", "#CD2626")
names(color_pallete) <- c("A", "C", "G", "T")

# removes outliers from vector
removeOutliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

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
    command <- paste(ClustalO, '-i', infile, '-o', outfile,'--force' ,threads, format, iterations)
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

# pattern matching in sequences, plots sequence logo around pattern 
signalLogo <- function(sequences, signal_pattern, LTR_sample_class){
  
  # get the most common signal start 
  signal_start <- 
    signal_start_median %>% 
    filter(LTR_class == LTR_sample_class) %$% 
    start_med
  
  # get sequences end (either signal start + 20 or sequence end)
  sub_sequences_end <- ifelse(width(sequences) < signal_start + seq_width, width(sequences), signal_start + seq_width)
  
  # get the subsequence 40 nt around the signal start
  sub_sequences <- subseq(x = sequences, start = signal_start - seq_width, end = sub_sequences_end)
  
  # write .fasta 
  write.fasta(as.list(sub_sequences),
              nbchar = 80,
              names = names(sub_sequences),
              as.string = TRUE,
              file.out = paste0(file_name, "_signal.fasta"),
              open = "w")
  
  # align sequences by Clustal Omega
  sub_sequences_aligned <- ClustalO(infile = paste0(getwd(), "/", file_name, "_signal.fasta"), 
                                    outfile = paste0(getwd(), "/", file_name, "_signal_msa.fasta"), 
                                    threads = 8, 
                                    force = T, 
                                    iterations = logo_alignment_iterations)
  
  # plot alignment
  msaPlot(aligned_sequences = sub_sequences_aligned, 
          plot_file_name = paste0(file_name, "_signal_full_msa_plot.pdf"))
  
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
    
    # plot alignment with no gaps
    msaPlot(aligned_sequences = sub_sequences_aligned, 
            plot_file_name = paste0(file_name, "_signal_full_noGaps_msa_plot.pdf"))
    
  }
  
  # get the most common signal start in aligned data 
  signal_start <- 
    as.data.frame(vmatchPattern(signal_pattern, sub_sequences_aligned, max.mismatch = pattern_mismatch)) %>%
    dplyr::count(start) %>% 
    arrange(desc(n)) %>% 
    dplyr::filter(start > 20) %>% 
    .[1, 1] %>% 
    as.integer()
  
  if(length(signal_start) > 0){
    # write logo with fixed signal start 
    weblogo(seqs = as.character(sub_sequences_aligned),
            open = F, 
            file.out = paste0(getwd(), "/", file_name, "_signal_msa_logo.pdf"),
            format = "pdf", 
            color.scheme = "classic",
            lower = signal_start - 20, 
            upper = signal_start + 5, 
            show.xaxis = F, 
            errorbars = T)
    
    # plot MSA of signal
    sub_sequences <- subseq(x = sub_sequences_aligned, start = signal_start - 20, end = signal_start + 5)
    msaPlot(aligned_sequences = sub_sequences, 
            plot_file_name = paste0(file_name, "_signal_short_msa_plot.pdf"))
    
    # print message
    cat(paste0(file_name, " logo done\n\n"))  
  } else cat(paste0(file_name, " logo not done\n\n"))
  
}

################################################################################## reading data
# LTR data
LTR_data <- 
  read_csv("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/seq_logo/LTR_ERVL_data_table.csv") %>% 
  dplyr::select(LTR_class, element_name = RepeatMasker_ID, consensus_width = consensus_sequence_length, consensus_length_fragment_numbers)

# repeatMasker
rptmsk <- read_delim("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_repeatMasker_mm10_20161012.txt.gz", delim =  "\t")

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

################################################################################## 
# get median of AATAAA pattern start in all LTR classes
signal_start_median <- 
  rptmsk_filter %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T) %>% 
  getSeq(x = Mmusculus, .) %>% 
  setNames(rptmsk_filter$fullName) %>% 
  vmatchPattern("AATAAA", .) %>% 
  unlist() %>% 
  as.data.frame() %>%
  separate(names, c("LTR_range", "LTR_subclass", "LTR_class"), sep = "\\|") %>% 
  group_by(LTR_class) %>% 
  summarize(start_med = round(median(start))) %>% 
  mutate(start_med = replace(start_med, LTR_class == "MTE", 170))

################################################################################## 
# for each LTR class:
# - get 200 random full length (consensus length +- 5%) LTR sequences
# - align them with ClustalO and plot as "heatmap"
# - find polyA signal and plot as sequence logo (with 20 nt downstream)

# parameters
seq_width <- 20
pattern_mismatch <- 1
logo_alignment_iterations <- 100

for(LTR_sample in unique(LTR_data$LTR_class)){
  
  # get 200 random LTRs from one class
  set.seed(1234)
  LTR_rptmsk_sample <- 
    rptmsk_filter %>% 
    dplyr::filter(LTR_class == LTR_sample) %>% 
    sample_n(ifelse(nrow(.) > 200, 200, nrow(.))) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
  
  # set output file names
  file_name <- paste0(LTR_sample, "_random_", length(LTR_rptmsk_sample))
  
  # get sequences 
  LTR_seqences <- getSeq(x = Mmusculus, LTR_rptmsk_sample)
  names(LTR_seqences) <- LTR_rptmsk_sample$fullName
  
  # add width, mismatch and alignment iteration to file names
  file_name <- paste0(file_name, 
                      "_width", seq_width, 
                      "_miss", pattern_mismatch, 
                      "_iter", logo_alignment_iterations)
  
  # get signal logo and MSA plot
  signalLogo(sequences = LTR_seqences, signal_pattern = "AATAAA", LTR_sample_class = LTR_sample)
  
}




