### INFO: 
### DATE: Fri Nov 08 21:22:18 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/PhD/algorithms_and_programming/2019_10_28/hackaton")

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
library(msa)
library(seqinr)
library(ape)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# subseq DNA string by start position
subseqStart <- function(dna_seq, position){
  
  subseq(dna_seq, start = position)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- file.path(getwd(), "prepare_data")

# set outpath
outpath <- getwd()

# animals table path
animals_path <- file.path(inpath, "animal_list.csv")

# gene coordinates path
genes_tb_path <- file.path(inpath, "hemoglobin.genomic_coordinates.all.csv")

# exonic coordinates path
exons_tb_path <- file.path(inpath, "hemoglobin.exonic_coordinates.all.csv")

# gene sequences path
genes_seq_path <- file.path(inpath, "hemoglobin.genomic_sequences.DNAStringSet.all.RDS")

# mouse hemoglobin path
mouse_hemoglobin_path <- file.path(inpath, "D1MQ2.hemoglobin.hsapiens.fasta")

######################################################## READ DATA
# read info about animals
animals_tb <- readr::read_csv(animals_path)

# read gene coordinates
genes_tb <- readr::read_csv(genes_tb_path)

# read exonic coordinates path
exons_tb <- readr::read_csv(exons_tb_path) 

# read gene sequences 
genes_seq <- readRDS(genes_seq_path)
 
# read mouse hemoglobin protein sequence
mouse_hemoglobin <- Biostrings::readAAStringSet(mouse_hemoglobin_path)

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
                seqnames = str_c(ensembl_name, gene_id, sep = ".")) %>% 
  dplyr::select(seqnames, start = exon_start, end = exon_end) %>% 
  dplyr::arrange(seqnames, start) %>% 
  dplyr::group_by(seqnames) %>% 
  dplyr::mutate(exon_id = str_c(seqnames, 1:n(), sep = ".")) %>%
  GRanges(.)

# extract exonic sequences
exons_seq <- genes_seq[exons_gr]
names(exons_seq) <- mcols(exons_gr)$exon_id

# join exons to spliced transcripts
transcripts_tb <- 
  tibble(exon_id = names(exons_seq), 
         exon_seq = as.character(exons_seq)) %>% 
  dplyr::mutate(transcript_id = str_remove(exon_id, "\\.\\d")) %>%
  dplyr::group_by(transcript_id) %>% 
  dplyr::summarise(transcript_seq = str_c(exon_seq, collapse = "")) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(ensembl_name = str_remove(transcript_id, "\\..*")) %>% 
  dplyr::left_join(., animals_tb %>% dplyr::select(ensembl_name, common_name), by = "ensembl_name") %>% 
  dplyr::mutate(transcript_name = str_c(common_name, str_remove(transcript_id, str_c(ensembl_name, ".")), sep = " "))


### find reading frames
# create an DNAStringSet
transcripts_seq <- DNAStringSet(transcripts_tb$transcript_seq)
names(transcripts_seq) <- transcripts_tb$transcript_name

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
  seq_protein <- Biostrings::translate(transcripts_all_rf[[seq_name]], if.fuzzy.codon = "solve")

  # get the longest protein
  longest_protein <- 
    seq_protein %>% 
    as.character(.) %>% 
    str_split(., "\\*") %>% 
    unlist(.) %>% 
    .[nchar(.) > 50] %>% 
    AAStringSet(.)
  
  # set names
  names(longest_protein) <- str_c(seq_name, 1:length(longest_protein), sep = ".")
  
  # return
  return(longest_protein)
  
}) %>% 
  set_names(., names(transcripts_all_rf))


### align predicted proteins with mouse hemoglobin subunit alpha
# iterate over sequences
hemoglobin_seq <- lapply(names(proteins_seq), function(seq_name){
  
  cat(seq_name, "\n")
  
  # extract sequence
  protein_seq <- proteins_seq[[seq_name]]
  
  # do the pairwise alignment with mouse hemoglobin
  alignment_scores <- sapply(names(protein_seq), function(name){
    
    # align and get score, 
    alignment_score <- 
      Biostrings::pairwiseAlignment(pattern = mouse_hemoglobin, subject = protein_seq[name]) %>% 
      score(.)
    
    # return
    return(alignment_score)
    
  }) %>% 
    set_names(names(protein_seq))
  
  # get protein with highest score
  best_protein <- protein_seq[which.max(alignment_scores)]

  # add score to best protein
  names(best_protein) <- str_c(names(best_protein), round(alignment_scores[names(best_protein)], 2), sep = "|")
  
  # return
  return(best_protein)
  
}) %>% 
  do.call(c, .)


### multiple sequence alignment and phylogeny tree
# MSA2
hemoglobin_msa <- msa::msa(hemoglobin_seq)

# convert to seqinr format
hemoglobin_msa_seqinr <- msaConvert(hemoglobin_msa, type = "seqinr::alignment")

# calculate distance matrix
hemoglobin_msa_dist <- seqinr::dist.alignment(hemoglobin_msa_seqinr, "identity")

# get tree and plot
hemoglobin_msa_tree <- ape::nj(hemoglobin_msa_dist)
plot(hemoglobin_msa_tree)


