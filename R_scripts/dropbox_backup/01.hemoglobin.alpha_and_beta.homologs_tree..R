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
inpath <- getwd()

# set outpath
outpath <- getwd()

# animals table path
hemoglobin_path <- file.path(inpath, "hemoglobin.alpha_and_beta.homologs_ensembl.csv")

# gene coordinates path
genes_tb_path <- file.path(inpath, "hemoglobin.genomic_coordinates.csv")

# exonic coordinates path
exons_tb_path <- file.path(inpath, "hemoglobin.exonic_coordinates.csv")

# gene sequences path
genes_seq_path <- file.path(inpath, "hemoglobin.genomic_sequences.DNAStringSet.RDS")

# human hemoglobin sequences
human_hemoglobin_path <- file.path(inpath, "hemoglobin.human_sequences.fasta")

######################################################## READ DATA
# read info about animals
hemoglobin_tb <- readr::read_csv(hemoglobin_path)

# read gene coordinates
genes_tb <- readr::read_csv(genes_tb_path)

# read exonic coordinates path
exons_tb <- readr::read_csv(exons_tb_path) 

# read gene sequences 
genes_seq <- readRDS(genes_seq_path)

# read human hemoglobin protein sequence
human_hemoglobin_seq <- Biostrings::readAAStringSet(human_hemoglobin_path)

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
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., hemoglobin_tb, by = "gene_id") %>% 
  dplyr::mutate(transcript_name = str_c(common_name, "hemoglobin", hemoglobin, sep = " "))


### find reading frames
# create an DNAStringSet
transcripts_seq <- DNAStringSet(transcripts_tb$transcript_seq)
names(transcripts_seq) <- transcripts_tb$transcript_name


library(systemPipeR)

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


### align predicted proteins with human hemoglobin subunit alpha/beta
# iterate over sequences
hemoglobin_seq <- lapply(names(proteins_seq), function(seq_name){
  
  cat(seq_name, "\n")
  
  # extract sequence
  protein_seq <- proteins_seq[[seq_name]]
  
  # extract human sequence human sequence
  human_seq <- human_hemoglobin_seq[str_detect(names(human_hemoglobin_seq), str_extract(seq_name, "alpha|beta"))]
  
  # get score
  alignment_scores <- Biostrings::pairwiseAlignment(pattern = protein_seq, 
                                                    subject = rep(human_seq, length(protein_seq)), 
                                                    scoreOnly = T)
  
  # get protein with highest score
  best_protein <- protein_seq[which.max(alignment_scores)]

  # return
  return(best_protein)
  
}) %>% 
  do.call(c, .)


### multiple sequence alignment and phylogeny tree
# MSA
hemoglobin_msa <- msa::msa(hemoglobin_seq)

# convert to seqinr format
hemoglobin_msa_seqinr <- msaConvert(hemoglobin_msa, type = "seqinr::alignment")

# calculate distance matrix
hemoglobin_msa_dist <- seqinr::dist.alignment(hemoglobin_msa_seqinr, "identity")

# get tree and plot
hemoglobin_msa_tree <- ape::nj(hemoglobin_msa_dist)
hemoglobin_msa_tree <- ape::ladderize(hemoglobin_msa_tree)
plot(hemoglobin_msa_tree)





#####################
### multiple sequence alignment and phylogeny tree - genomic
# add names to genomic sequences
genes_seq_tb <- 
  tibble(gene_seq = as.character(genes_seq), 
         gene_id = names(genes_seq)) %>% 
  dplyr::left_join(., hemoglobin_tb, by = "gene_id") %>% 
  dplyr::mutate(gene_name = str_c(common_name, hemoglobin, sep = " "))

# named sequences
genes_seq_named <- DNAStringSet(genes_seq_tb$gene_seq)
names(genes_seq_named) <- genes_seq_tb$gene_name

# MSA
hemoglobin_msa <- msa::msa(genes_seq_named)

# convert to seqinr format
hemoglobin_msa_seqinr <- msaConvert(hemoglobin_msa, type = "seqinr::alignment")

# calculate distance matrix
hemoglobin_msa_dist <- seqinr::dist.alignment(hemoglobin_msa_seqinr, "identity")

# get tree and plot
hemoglobin_msa_tree <- ape::njs(hemoglobin_msa_dist)
hemoglobin_msa_tree <- ape::ladderize(hemoglobin_msa_tree)
plot(hemoglobin_msa_tree)


### multiple sequence alignment and phylogeny tree - transcriptomic
# add names to genomic sequences
# named sequences
transcipt_seq_named <- DNAStringSet(transcripts_tb$transcript_seq)
names(transcipt_seq_named) <- transcripts_tb$transcript_name

# MSA
hemoglobin_msa <- msa::msa(transcipt_seq_named)

# convert to seqinr format
hemoglobin_msa_seqinr <- msaConvert(hemoglobin_msa, type = "seqinr::alignment")

# calculate distance matrix
hemoglobin_msa_dist <- seqinr::dist.alignment(hemoglobin_msa_seqinr, "identity")

# get tree and plot
hemoglobin_msa_tree <- ape::njs(hemoglobin_msa_dist)
hemoglobin_msa_tree <- ape::ladderize(hemoglobin_msa_tree)
plot(hemoglobin_msa_tree)

predORF
