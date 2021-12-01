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

changeWD("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/LTR_alignment/seq_logo/exapted_LTRs")

################################################################################## functions
# mannualy set colors to the stages
# (A) green; (C) blue; (G) yellow; (T) red
color_pallete <- c("#66CD00", "#104E8B", "#FFD700", "#CD2626")
names(color_pallete) <- c("A", "C", "G", "T")

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
    command <- paste(ClustalO, '-i', infile, '-o', outfile, '--force', threads, format, iterations)
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
          #           axis.text.x = element_text(size = 3, hjust = 1, vjust = 0, angle = 90),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.title.y = element_blank(),
          #           axis.text.y = element_text(size = 1, hjust = 1, vjust = 0),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggsave(plot_file_name)
  
  cat(paste0(plot_file_name, " done\n"))
  
}

################################################################################## reading data
# Vedran's figure 3 data 
repeats_list <- readRDS("/common/WORK/vfranke/Projects/PSvoboda_MT/Results/MT_FindRepeats/mm/160907.mm.Results.rds", refhook = NULL)

# LTR order 
LTR_order <- c("MLT1", "MLT2", 
               "ORR1G", "ORR1F", "ORR1E", "ORR1D", "ORR1C", "ORR1B", "ORR1A", 
               "MTE", "MTE2", "MTD", "MTC", "MTB", "MTA", "MT2A", "MT2B", "MT2C", "MT2")

# LTR data with consensus length (Maja)
consensus_length <- 
  read_csv("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/LTR_alignment/seq_logo/LTR_ERVL_data_table.csv") %>% 
  dplyr::select(LTR_class, LTR_subclass = RepeatMasker_ID, consensus_width = consensus_sequence_length) %>% 
  mutate(LTR_class = replace(LTR_class, LTR_class == "MTB2", "MTB")) %>% 
  group_by(LTR_class) %>%
  summarize(consensus_width = max(consensus_width)) %>% 
  mutate(LTR_class = factor(LTR_class, levels = LTR_order)) 

################################################################################## read annotated exons' splice donor sites
# read ENSEMBL and UCSC splice donors
splice_donors <- 
  rbind(read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/splicing/Ensembl_GRCm38.86.20161128_splice_donor.txt", delim =  "\t"), 
        read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/splicing/UCSC_knownGene_mm10_20161126_splice_donor.txt", delim =  "\t")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) 

################################################################################## read splice junctions coordinates
# get splice junction starts in GV
all_strands <- setNames(c("*", "+", "-"), as.character(0:2))
all_intron_motifs <- setNames(c("non-canonical", "GT/AG", "CT/AC", "GC/AG", "CT/GC", "AT/AC", "GT/AT"), as.character(0:6))

sj_splice_start <- 
  read_delim("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_GV.WE/s_GV.WESJ.out.tab", 
             delim = "\t", col_names = c("seqnames", "start", "end", "strand", "intron_motif", "annot", "n_uniq", "n_multi", "overhang")) %>%
  mutate(strand = all_strands[as.character(strand)], 
         intron_motif = all_intron_motifs[as.character(intron_motif)], 
         SJ_fullName = paste0(seqnames, ":",
                              start, "-", 
                              end, ":", 
                              strand)) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T) %>% 
  GenomicRanges::resize(., fix = "start", width = 1)

################################################################################## 
# full and partial contribution to a 5’ exon in protein-coding genes and lncRNA
selected_LTRs <- 
  repeats_list$RepeatsSelected_MALR %>% 
  filter(ex.category == "5' exon", 
         !is.na(repClassCust)) %>% 
  dplyr::select(coordinates = rep.coord, strand = rep.strand, LTR_subclass = repName, LTR_class = repClassCust) %>% 
  left_join(consensus_length, by = "LTR_class") %>% 
  mutate(fullName = paste0(coordinates, ":", 
                           strand, "|", 
                           LTR_subclass, "|",
                           LTR_class)) %>% 
  separate(coordinates, c("seqnames", "coordinates"), ":") %>% 
  separate(coordinates, c("start", "end"), "-") %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

################################################################################## get splice start distribution
# overlap splice junction starts with selected LTRs
sj_splice_start_LTR <- sj_splice_start[subjectHits(findOverlaps(selected_LTRs, sj_splice_start))] 
sj_splice_start_LTR$fullName <- selected_LTRs[queryHits(findOverlaps(selected_LTRs, sj_splice_start))]$fullName

# # filter out MTA/MTB
# sj_splice_start_LTR <- sj_splice_start_LTR[!grepl(pattern = "MTA|MTB", x = sj_splice_start_LTR$fullName)]

# overlap splice starts with annotated splice donors 
exon_splice_donors <- sj_splice_start_LTR[queryHits(findOverlaps(sj_splice_start_LTR, splice_donors))]$SJ_fullName

# filter out those overlaping with annotated exons, expand ranges to +-10 nt
expand_range <- 10
sj_splice_start_LTR %<>%  
  as.data.frame() %>% 
  dplyr::filter(!(SJ_fullName %in% exon_splice_donors)) %>% 
  dplyr::select(c(1:3, 5:6, 12)) %>% 
  mutate(start = start - expand_range, 
         end = end + expand_range) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

# get splice site sequences
LTR_seqences <- 
  getSeq(x = Mmusculus, sj_splice_start_LTR) %>% 
  setNames(sj_splice_start_LTR$fullName)

# set file names
file_name <- paste0("spliceSite_", expand_range, "_MT2")

# # write .fasta 
# write.fasta(as.list(LTR_seqences),
#             nbchar = 80,
#             names = names(LTR_seqences),
#             as.string = TRUE,
#             file.out = paste0(file_name, ".fasta"),
#             open = "w")

# align sequences by Clustal Omega
sequences_aligned <- ClustalO(infile = paste0(getwd(), "/", file_name, ".fasta"), 
                              outfile = paste0(getwd(), "/", file_name, "_msa.fasta"), 
                              threads = 8, 
                              force = F, 
                              iterations = 100)

# # plot alignment
# msaPlot(aligned_sequences = sequences_aligned, 
#         plot_file_name = paste0(file_name, "_msa_plot.pdf"))


# get gap index
gap_index <- 
  consensusMatrix(sequences_aligned, as.prob = F, baseOnly = T)[5, ] %>%
  is_greater_than(length(sequences_aligned) * 0.8) %>% 
  which()

if(length(gap_index) > 0){
  
  # keep only those columns in alignment
  sequences_aligned <- 
    as.matrix(sequences_aligned) %>% 
    .[, -gap_index] %>% 
    apply(., 1, paste, collapse = "") %>% 
    DNAStringSet()
  
  # # plot alignment with no gaps
  # msaPlot(aligned_sequences = sequences_aligned, 
  #         plot_file_name = paste0(file_name, "_noGaps_msa_plot.pdf"))
  # 
}

# write logo with fixed signal start 
weblogo(seqs = as.character(sequences_aligned),
        open = F, 
        file.out = paste0(getwd(), "/test/", file_name, "_signal_msa_logo.pdf"),
        format = "pdf", 
        color.scheme = "classic",
        show.xaxis = F, 
        errorbars = F)

