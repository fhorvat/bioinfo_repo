rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd(".")

######################################################## LIBRARIES
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(reshape2)
library(magrittr)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(BiocParallel)
library(BSgenome.Mmusculus.UCSC.mm10)
library(seqinr)
library(Biostrings)
library(RWebLogo)

######################################################## FUNCTIONS
# Clustal Omega function
ClustalO <- function(infile, outfile, ClustalO = NULL, threads = 12, 
                     what = 'DNA', force = T, format = 'fa', iterations = 0){
  
  if(!what %in% c('AA', 'DNA'))
    stop('can only align DNA or AA')
  if(is.null(ClustalO))
    ClustalO <- 'clustal-omega-1.2.3/bin/clustalo'
  
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

######################################################## READ DATA
# LTR order 
LTR_order <- c("MLT1", "MLT2", 
               "ORR1G", "ORR1F", "ORR1E", "ORR1D", "ORR1C", "ORR1B", "ORR1A", 
               "MTE", "MTE2", "MTD", "MTC", "MTB", "MTA", "MT2A", "MT2B", "MT2C", "MT2")

# LTR data consensus length (2016_Franke_Supplemental_table_S1)
consensus_length <- 
  read_csv("LTR_ERVL_data_table.csv") %>% 
  dplyr::select(LTR_class, LTR_subclass = RepeatMasker_ID, consensus_width = consensus_sequence_length) %>% 
  mutate(LTR_class = replace(LTR_class, LTR_class == "MTB2", "MTB")) %>% 
  group_by(LTR_class) %>%
  summarize(consensus_width = max(consensus_width)) %>% 
  mutate(LTR_class = factor(LTR_class, levels = LTR_order)) 

# ENSEMBL and UCSC splice donors coordinates 
splice_donors <- 
  rbind(read_delim("Ensembl_GRCm38.86.20161128_splice_donor.txt", delim =  "\t"), 
        read_delim("UCSC_knownGene_mm10_20161126_splice_donor.txt", delim =  "\t")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) 

# splice junction starts in GV stage (SJ.out.tab from STAR RNA-seq aligner)
all_strands <- setNames(c("*", "+", "-"), as.character(0:2))
all_intron_motifs <- setNames(c("non-canonical", "GT/AG", "CT/AC", "GC/AG", "CT/GC", "AT/AC", "GT/AT"), as.character(0:6))

sj_splice_start <- 
  read_delim("s_GV.WE.SJ.out.tab", 
             delim = "\t", col_names = c("seqnames", "start", "end", "strand", "intron_motif", "annot", "n_uniq", "n_multi", "overhang")) %>%
  mutate(strand = all_strands[as.character(strand)], 
         intron_motif = all_intron_motifs[as.character(intron_motif)], 
         SJ_fullName = paste0(seqnames, ":",
                              start, "-", 
                              end, ":", 
                              strand)) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T) %>% 
  GenomicRanges::resize(., fix = "start", width = 1)

# LTRs with full and partial contribution to a 5’ exon in protein-coding genes and lncRNA (figure 3 data)
repeats_list <- readRDS("160907.mm.Results.rds", refhook = NULL)

all_LTRs <- 
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
  dplyr::mutate(start = as.integer(start), 
                end = as.integer(end))

######################################################## MAIN CODE
### splice site motif
LTR_families <- c("all", "MT", "MT2", "ORR")

invisible(lapply(LTR_families, function(x){
  
  # set file names
  file_name <- paste0("spliceSite_", x)
  
  # select LTRs by family
  if(x == "all"){
    selected_LTRs <- 
      all_LTRs %>% 
      dplyr::filter(!(LTR_class %in% c("MTA", "MTB")))
  }else{
    if(x == "MT"){
      selected_LTRs <- 
        all_LTRs %>% 
        dplyr::filter(str_detect(LTR_class, "MT"), 
                      !str_detect(LTR_class, "MT2"))
    }else{
      selected_LTRs <- 
        all_LTRs %>% 
        dplyr::filter(str_detect(LTR_class, x))
    }
  }
  
  selected_LTRs <- makeGRangesFromDataFrame(selected_LTRs, keep.extra.columns = T)
  
  # overlap splice starts with selected LTRs
  sj_splice_start_LTR <- sj_splice_start[subjectHits(findOverlaps(selected_LTRs, sj_splice_start))] 
  sj_splice_start_LTR$fullName <- selected_LTRs[queryHits(findOverlaps(selected_LTRs, sj_splice_start))]$fullName
  
  # overlap splice starts with annotated splice donors 
  exon_splice_donors <- sj_splice_start_LTR[queryHits(findOverlaps(sj_splice_start_LTR, splice_donors))]$SJ_fullName
  
  # remove splice starts overlaping annotated splice donors
  sj_splice_start_LTR %<>%  
    as.data.frame() %>% 
    dplyr::filter(!(SJ_fullName %in% exon_splice_donors)) %>% 
    dplyr::select(c(1:3, 5:6, 12)) %>% 
    mutate(start = start - 10, 
           end = end + 10) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
  
  # splice site sequences
  LTR_seqences <- 
    getSeq(x = Mmusculus, sj_splice_start_LTR) %>% 
    setNames(sj_splice_start_LTR$fullName)
  
  # write .fasta
  write.fasta(as.list(LTR_seqences),
              nbchar = 80,
              names = names(LTR_seqences),
              as.string = TRUE,
              file.out = paste0(file_name, ".fasta"),
              open = "w")
  
  # align sequences with Clustal Omega
  sequences_aligned <- ClustalO(infile = paste0(file_name, ".fasta"), 
                                outfile = paste0(file_name, "_msa.fasta"), 
                                threads = 8, 
                                force = F, 
                                iterations = 100)
  
  # gap index
  gap_index <- 
    consensusMatrix(sequences_aligned, as.prob = F, baseOnly = T)[5, ] %>%
    is_greater_than(length(sequences_aligned) * 0.8) %>% 
    which()
  
  if(length(gap_index) > 0){
    sequences_aligned <- 
      as.matrix(sequences_aligned) %>% 
      .[, -gap_index] %>% 
      apply(., 1, paste, collapse = "") %>% 
      DNAStringSet()
  }
  
  # save logo
  weblogo(seqs = as.character(sequences_aligned),
          open = F, 
          file.out = paste0(file_name, "_splice_site_logo.pdf"),
          format = "pdf", 
          color.scheme = "classic",
          show.xaxis = F, 
          errorbars = F)
  
}))
