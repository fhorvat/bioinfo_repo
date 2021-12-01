### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/annotation/blast_consensus_L1")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS
expandRange = function(x, upstream=0, downstream=0) {
  strand_is_minus = strand(x) == "-"
  on_plus = which(!strand_is_minus)
  on_minus = which(strand_is_minus)
  start(x)[on_plus] = start(x)[on_plus] - upstream
  start(x)[on_minus] = start(x)[on_minus] - downstream
  end(x)[on_plus] = end(x)[on_plus] + downstream
  end(x)[on_minus] = end(x)[on_minus] + upstream
  x
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# rmsk path
rmsk_path <- list.files(path = genome_dir, pattern = "rmsk.*\\.joined_rmsk_id\\.fa\\.out\\.gz", full.names = T)

# blast results path
blast_results_path <- list.files(inpath, ".*\\.blastn.txt", full.names = T)

# consensus fasta path
consensus_fasta_path <- list.files(inpath, ".*\\.consensus\\.fasta", full.names = T)

######################################################## READ DATA
# read blast results
blast_results <- purrr::map(blast_results_path, function(path){
  
  # read table and name it
  readr::read_delim(file = path, delim = "\t", col_names = c("query_id", "subject_id", "identity_perc", "alignment_length", 
                                                             "mismatches", "gap_open", 
                                                             "query_start", "query_end", "subject_start", "subject_end", 
                                                             "e_value", "bit_score")) %>% 
    arrange(-alignment_length) 
  
}) %>% 
  magrittr::set_names(., basename(blast_results_path) %>% str_remove(., "\\.blastn\\.txt"))

# read consensus fasta
consensus_fasta <- 
  purrr::map(consensus_fasta_path, Biostrings::readDNAStringSet) %>% 
  do.call(c, .)
names(consensus_fasta) <- basename(consensus_fasta_path) %>% str_remove(., "\\.fasta")

# read repeatMasker
rmsk_gr <- 
  readr::read_delim(rmsk_path, delim = "\t") %>% 
  dplyr::filter(repClass == "LINE") %>% 
  dplyr::mutate(strand = "*") %>% 
  GRanges(.)

######################################################## MAIN CODE
l1_name <- "L1_full_length_manual_200707.without_5p_repeat.consensus"

# filter results
line1_tb <- 
  blast_results[[l1_name]] %>% 
  dplyr::mutate(strand = ifelse(subject_start < subject_end, "+", "-")) %>% 
  dplyr::mutate(start = ifelse(strand == "+", subject_start, subject_end), 
                end = ifelse(strand == "+", subject_end, subject_start)) %>% 
  dplyr::select(seqnames = subject_id, start, end, strand, query_start, query_end, identity_perc) %>% 
  dplyr::mutate(query_alignment_width = query_end - query_start + 1, 
                subject_alignment_width = end - start + 1, 
                rmsk_id = str_c(seqnames, ":", start, "-", end)) %>% 
  dplyr::filter(query_start == 1) %>% 
  dplyr::filter(query_alignment_width >= width(consensus_fasta[l1_name]) - 20)

# create GRanges object 
line1_gr <- 
  line1_tb %>% 
  GRanges(.)

# extract sequences
line1_seq <- Biostrings::getSeq(x = BSgenome.Maur.UCSC.Siomi, names = line1_gr)
names(line1_seq) <- line1_gr$rmsk_id

# find sequences containing Ns
n_sequences <- 
  vmatchPattern(pattern = "N", subject = line1_seq) %>% 
  unlist(.) %>% 
  names(.) %>% 
  unique(.)

# remove sequences containg N
line1_seq <- line1_seq[!names(line1_seq) %in% n_sequences]


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


# add info about ORFs to original table, filter
line1_tb_filt <- 
  line1_tb %>% 
  dplyr::left_join(., insertions_orfs_tb, by = "rmsk_id") %>% 
  dplyr::filter(rmsk_id %in% names(line1_seq)) %>% 
  dplyr::filter(longest_orf_1 >= 1200*3, 
                longest_orf_2 >= 370*3) %>% 
  dplyr::filter(longest_orf_1 == 3816, 
                longest_orf_2 == 1491)

# save
readr::write_csv(line1_tb_filt %>% dplyr::select(-c(query_alignment_width)), file.path(outpath, str_c(l1_name, "BLAST_hits", "with_both_ORFs", "csv", sep = ".")))


### expand ranges
# create GRanges, expand them 2000 upstream, get sequences, save
line1_filt_gr <- 
  line1_tb_filt %>% 
  GRanges(.)

# for each hit find rmsk_ID to which it belongs 
overlaps <- findOverlaps(line1_filt_gr, rmsk_gr, minoverlap = 5000)
hits_tb <- 
  rmsk_gr[subjectHits(overlaps)] %>% 
  as_tibble(.) %>% 
  dplyr::mutate(hit_coordinates = mcols(line1_filt_gr[queryHits(overlaps)])$rmsk_id) %>% 
  dplyr::group_by(hit_coordinates) %>%
  dplyr::filter(width == min(width)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::select(hit_coordinates, rmsk_id)

# add rmsk to table
line1_filt_gr %<>% 
  as_tibble(.) %>% 
  dplyr::rename(hit_coordinates = rmsk_id) %>% 
  dplyr::left_join(., hits_tb, by = "hit_coordinates") %>% 
  GRanges(.)

# expand
seqlevels(line1_filt_gr) <- seqlevels(BSgenome.Maur.UCSC.Siomi)
seqlengths(line1_filt_gr) <- seqlengths(BSgenome.Maur.UCSC.Siomi)
line1_gr_expand <- expandRange(line1_filt_gr, upstream = 2000, downstream = 0)
line1_gr_expand <- trim(line1_gr_expand)
names(line1_gr_expand) <- str_c(line1_gr_expand$hit_coordinates, line1_gr_expand$rmsk_id, sep = ".")

# save as bed
rtracklayer::export.bed(object = line1_gr_expand, con = file.path(outpath, str_c(l1_name, "BLAST_hits", "with_both_ORFs", "2kb_upstream", "bed", sep = ".")))


# get sequences
line1_seq <- getSeq(BSgenome.Maur.UCSC.Siomi, line1_gr_expand)
line1_seq <- line1_seq[!str_detect(line1_seq, "N")]

# save sequences
Biostrings::writeXStringSet(line1_seq, file = file.path(outpath, str_c(l1_name, "BLAST_hits", "with_both_ORFs", "2kb_upstream", "fasta", sep = ".")))
