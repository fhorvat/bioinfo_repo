### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/annotation/04_L1_ORFs_blast_hits")

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
library(BSgenome.Maur.UCSC.Siomi)
library(systemPipeR)
library(msa)
library(ape)
library(seqinr)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# rmsk paths
rmsk_joined_path <- list.files(path = genome_dir, pattern = "rmsk.*\\.joined_rmsk_id\\.fa\\.out\\.gz", full.names = T)
rmsk_clean_path <- list.files(path = genome_dir, pattern = "rmsk.*\\.clean\\.fa\\.out\\.gz", full.names = T)

# blast results path
blast_results_path <- list.files(inpath, ".*\\.blastn.txt", full.names = T)

# consensus fasta path
orfs_fasta_path <- file.path(inpath, "L1_full_length_manual_200707.consensus.ORFs.fasta")

######################################################## READ DATA
# read blast results
blast_results <- readr::read_delim(file = blast_results_path, 
                                   delim = "\t", 
                                   col_names = c("query_id", "subject_id", "identity_perc", "alignment_length", 
                                                 "mismatches", "gap_open", 
                                                 "query_start", "query_end", "subject_start", "subject_end", 
                                                 "e_value", "bit_score")) 

# read ORFs fasta
orfs_fasta <- Biostrings::readDNAStringSet(orfs_fasta_path)

# read joined repeatMasker
rmsk_joined_gr <- 
  readr::read_delim(rmsk_joined_path, delim = "\t") %>% 
  dplyr::filter(repClass == "LINE") %>% 
  dplyr::mutate(strand = "*") %>% 
  GRanges(.)

# read clean repeatMasker
rmsk_clean_gr <- 
  readr::read_delim(rmsk_clean_path, delim = "\t") %>% 
  dplyr::filter(repClass == "LINE") %>% 
  GRanges(.)

######################################################## MAIN CODE
# filter results
line1_tb <- 
  blast_results %>% 
  dplyr::mutate(strand = ifelse(subject_start < subject_end, "+", "-")) %>% 
  dplyr::mutate(start = ifelse(strand == "+", subject_start, subject_end), 
                end = ifelse(strand == "+", subject_end, subject_start)) %>% 
  dplyr::select(seqnames = subject_id, start, end, strand, query_start, query_end, identity_perc) %>% 
  dplyr::mutate(query_alignment_width = query_end - query_start + 1, 
                subject_alignment_width = end - start + 1, 
                rmsk_id = str_c(seqnames, ":", start, "-", end)) %>% 
  dplyr::filter(query_start <= 10) %>%
  dplyr::filter(query_end >= width(orfs_fasta) - 10)

# create GRanges object 
line1_gr <- 
  line1_tb %>% 
  GRanges(.)

# for each hit find rmsk_ID to which it belongs - when it picks more than one take "within" and not "interupted"
overlaps <- findOverlaps(line1_gr, rmsk_clean_gr, ignore.strand = T)
hits_tb <- 
  rmsk_clean_gr[subjectHits(overlaps)] %>% 
  as_tibble(.) %>% 
  dplyr::mutate(hit_coordinates = mcols(line1_gr[queryHits(overlaps)])$rmsk_id) %>% 
  dplyr::group_by(hit_coordinates, rmsk_id) %>%
  dplyr::summarise(width = sum(width)) %>% 
  dplyr::group_by(hit_coordinates) %>% 
  dplyr::filter(width == max(width)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., rmsk_joined_gr %>% as_tibble(.) %>% dplyr::select(rmsk_id, repName), by = "rmsk_id") %>% 
  dplyr::select(hit_coordinates, rmsk_id, repName)

# add rmsk ID to table
line1_gr %<>% 
  as_tibble(.) %>% 
  dplyr::rename(hit_coordinates = rmsk_id) %>% 
  dplyr::left_join(., hits_tb, by = "hit_coordinates") %>% 
  GRanges(.)

# extract sequences
line1_seq <- Biostrings::getSeq(x = BSgenome.Maur.UCSC.Siomi, names = line1_gr)
names(line1_seq) <- str_c(line1_gr$repName, ".", line1_gr$rmsk_id)
names(line1_gr) <- str_c(line1_gr$repName, ".", line1_gr$rmsk_id)


### check if sequences actually have ORFs
# find ORFs
line1_orfs <- predORF(line1_seq, n = "all", type = "grl", mode = "orf", strand = "both", longest_disjoint = F)

# join to one GRanges object
insertions_orfs_list <- 
  line1_orfs %>% 
  unname(.) %>% 
  do.call(c, .)

# set repeatMasker id
mcols(insertions_orfs_list)$rmsk_id <- as.character(seqnames(insertions_orfs_list))

# split by name and reading frame
insertions_orfs_list <- split(insertions_orfs_list, list(insertions_orfs_list$rmsk_id, insertions_orfs_list$inframe2end))


### reduce ORF coordinates per reading frame and strand in each insertion
# loop
insertions_orfs_reduced <- purrr::map(names(insertions_orfs_list), function(orf_name){
  
  # get one insertion/frame
  insertions_orf <- insertions_orfs_list[[orf_name]]
  
  # get name and reading frame
  rmsk_id <- unique(mcols(insertions_orf)$rmsk_id)
  
  # to GRanges, reduce
  insertions_orf <- 
    insertions_orfs_list[[orf_name]] %>% 
    GenomicRanges::reduce(., ignore.strand = F, min.gapwidth = 0)
  
  # set rmsk_id and reading frames
  mcols(insertions_orf)$rmsk_id <- rmsk_id
  
  # return 
  return(insertions_orf)
  
})

### find two longest ORFs per rmsk_id and save
# get tidy table
insertions_orfs_tb <- 
  insertions_orfs_reduced %>% 
  do.call(c,.) %>% 
  as_tibble(.)

insertions_orfs_tb %<>% 
  dplyr::select(rmsk_id, start, end, strand, width) %>%
  dplyr::arrange(-width) %>% 
  dplyr::group_by(rmsk_id) %>% 
  dplyr::summarise(longest_orf_1 = sort(width, decreasing = T)[1], 
                   longest_orf_2 = sort(width, decreasing = T)[2]) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate_at(vars(dplyr::starts_with("longest_orf")), ~replace(., is.na(.), 0)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(longest_orf_1 >= 1200*3, 
                longest_orf_2 >= 370*3)

# get tidy table
line1_orfs_tb <- insertions_orfs_tb

# filter
line1_gr <- line1_gr[names(line1_gr) %in% line1_orfs_tb$rmsk_id]
line1_seq <- line1_seq[names(line1_seq) %in% line1_orfs_tb$rmsk_id]

# save table
line1_gr %>% 
  as_tibble(.) %>% 
  dplyr::select(-c(width, hit_coordinates)) %>%
  tidyr::unite(hit_coordinates, seqnames, start, end, sep = " ") %>% 
  readr::write_csv(file.path(outpath, str_c("L1_full_length_manual_200707.consensus.ORFs_blast_hits", "csv", sep = ".")))

# save sequences
Biostrings::writeXStringSet(line1_seq, file = file.path(outpath, str_c("L1_full_length_manual_200707.consensus.ORFs_blast_hits", "fasta", sep = ".")))


### MSA and phylogeny tree
# msa
seq_msa <- msa::msa(line1_seq, method = "ClustalOmega")

# write as fasta
seq_msa@unmasked %>%
  Biostrings::writeXStringSet(., file = file.path(outpath, str_c("L1_full_length_manual_200707.consensus.ORFs_blast_hits", "ClustalO.msa", "fasta", sep = ".")))

# convert to seqinr format
msa_seqinr <- msaConvert(seq_msa, type = "seqinr::alignment")

# calculate distance matrix
msa_dist <- seqinr::dist.alignment(msa_seqinr, "identity")

# get tree and plot
msa_tree <- ape::njs(msa_dist)
msa_tree <- ape::ladderize(msa_tree)

# save with labels
pdf(file = file.path(outpath, str_c("L1_full_length_manual_200707.consensus.ORFs_blast_hits", "ClustalO.msa", "phylo_tree", "pdf", sep = ".")), width = 20, height = 75)
plot.phylo(msa_tree, edge.width = 2, font = 1, show.tip.label = TRUE, cex = 1)
axisPhylo(backward = F)
dev.off()

