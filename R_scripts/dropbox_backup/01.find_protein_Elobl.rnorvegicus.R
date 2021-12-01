### INFO: 
### DATE: Fri Nov 08 21:22:18 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/Elob_MSA/rat_Elobl")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(purrr)

library(GenomicRanges)
library(Biostrings)
library(systemPipeR)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# subseq DNA string by start position
subseqStart <- function(dna_seq, position){
  
  subseq(dna_seq, start = position)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# gene coordinates path
genes_tb_path <- file.path(inpath, "Elobl.genomic_coordinates.csv")

# exonic coordinates path
exons_tb_path <- file.path(inpath, "Elobl.exonic_coordinates.csv")

# gene sequences path
genes_seq_path <- file.path(inpath, "Elobl.genomic_region.rnorvegicus.fa")

######################################################## READ DATA
# read gene coordinates
genes_tb <- readr::read_csv(genes_tb_path)

# read exonic coordinates path
exons_tb <- readr::read_csv(exons_tb_path) 

# read gene sequences 
genes_seq <- Biostrings::readDNAStringSet(genes_seq_path)
names(genes_seq) <- unique(genes_tb$gene_id)

######################################################## MAIN CODE
# prepare data
genes_tb %<>% 
  dplyr::rename(gene_start = start, gene_end = end)

### splice genes
# add gene coordinates to exons, calculate relative exon coordinates in genes
exons_gr <- 
  exons_tb %>% 
  dplyr::left_join(., genes_tb %>% dplyr::select(gene_id, gene_start, gene_end), by = "gene_id") %>% 
  dplyr::mutate(exon_start = start - gene_start + 1, 
                exon_end = end - gene_start + 1, 
                strand = "*", 
                seqnames = gene_id) %>% 
  dplyr::select(seqnames, start = exon_start, end = exon_end) %>% 
  dplyr::arrange(seqnames, start) %>% 
  dplyr::group_by(seqnames) %>% 
  dplyr::mutate(exon_id = str_c(seqnames, 1:n(), sep = ".")) %>%
  GRanges(.)

# # isto radi
# exons_seq <- getSeq(x = genes_seq, exons_gr)

# extract exonic sequences
exons_seq <- genes_seq[exons_gr]
names(exons_seq) <- mcols(exons_gr)$exon_id

# join exons to spliced transcripts
transcripts_tb <- 
  tibble(exon_id = names(exons_seq), 
         exon_seq = as.character(exons_seq)) %>% 
  dplyr::arrange(exon_id) %>% 
  dplyr::mutate(gene_id = str_remove(exon_id, "\\.\\d")) %>%
  dplyr::group_by(gene_id) %>% 
  dplyr::summarise(transcript_seq = str_c(exon_seq, collapse = "")) %>% 
  dplyr::ungroup(.)


### find reading frames
# create an DNAStringSet
transcripts_seq <- DNAStringSet(transcripts_tb$transcript_seq)
names(transcripts_seq) <- "Elobl.genome_browser"

# # save as fasta
# Biostrings::writeXStringSet(x = transcripts_seq, filepath = file.path(outpath, "Elobl.mRNA.rnorvegicus_predicted.fa"))

# read fasta from transcriptome assembly
trans_assembly <- Biostrings::readDNAStringSet(filepath = file.path(inpath, "blastn_transcriptome_assembly", "Elobl.mRNA.rnorvegicus.TRINITY_GG_981_c0_g1_i3.fa"))
names(trans_assembly) <- "Elobl.transcriptome_assembly"

# join with original sequence 
transcripts_seq <- c(transcripts_seq, trans_assembly)

# create a reverse complement
transcripts_rc <- reverseComplement(transcripts_seq)

# iterate over sequences to get all 6 reading frames for each sequence 
transcripts_all_rf <- lapply(names(transcripts_seq), function(seq_name){
  
  # iterate over 3 possible reading frames
  seq_rf_list <- lapply(1:3, function(frame){
    
    # subseq in 3 positions
    seq_frames <- subseqStart(transcripts_seq[[seq_name]], position = frame)
    seq_rc_frames <- subseqStart(transcripts_rc[[seq_name]], position = frame)
    
    # join to one list and return
    return(list(seq_frames, seq_rc_frames))
    
  }) %>% 
    unlist(.)
  
  # return a DNAStringSet list 
  return(DNAStringSet(seq_rf_list))
  
}) %>% 
  set_names(., names(transcripts_seq))


### translate to proteins, get predicted proteins longer than 100
# iterate over sequences
proteins_seq <- lapply(names(transcripts_all_rf), function(seq_name){
  
  cat(seq_name, "\n")
  
  # translate
  seq_protein <- Biostrings::translate(x = transcripts_all_rf[[seq_name]], 
                                       genetic.code = GENETIC_CODE, 
                                       if.fuzzy.codon = "solve", 
                                       no.init.codon = T)
  
  # get proteins starting with M
  longest_protein <- 
    seq_protein %>% 
    as.character(.) %>% 
    str_split(., "\\*") %>% 
    unlist(.) %>% 
    .[str_detect(., "M")] %>% 
    str_remove(., ".*?(?=M)") %>% 
    .[width(.) > 50]
  
  # if not empty
  if(length(longest_protein) > 0){
    
    # create AAStringSet
    longest_protein %<>% AAStringSet(.)
    
    # set names
    names(longest_protein) <- str_c(seq_name, 1:length(longest_protein), sep = ".")
    
  } 
  
  # return
  return(longest_protein)
  
}) %>% 
  set_names(., names(transcripts_all_rf))

# save as protein fasta
final_seq <- proteins_seq$Elobl.genome_browser
names(final_seq) <- "Elobl rnorvegicus transcriptome assembly"
Biostrings::writeXStringSet(final_seq, file.path(outpath, "../sequences/Elobl.protein.rnorvegicus.transcriptome_assembly.fa"))
