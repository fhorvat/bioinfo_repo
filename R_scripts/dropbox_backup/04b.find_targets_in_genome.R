### INFO: 
### DATE: Thu Jul 18 19:52:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/miR-205_pig/genomic_targets")

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

library(Biostrings)
library(GenomicRanges)
library(GenomicAlignments)
library(BSgenome.Sscrofa.UCSC.susScr11)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 100


# ### exon coordinates
# # genome path
# genome_dir <- "/common/DB/genome_reference/pig/susScr11.Sscrofa11.1.GCA_000003025.6"
# 
# # gene info path
# genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames\\.geneInfo\\.csv$"), full.names = T)


### FPKM
# FPKM table path
fpkm_tb_path <- list.files(path = file.path(inpath, "expression"), pattern = str_c(".*FPKM_mean\\.csv$"), full.names = T)

######################################################## READ DATA
# # read transcripts info
# genes_info <- readr::read_csv(genes_info_path)

# read FPKM table
fpkm_tb <- readr::read_csv(fpkm_tb_path)

######################################################## MAIN CODE
### get covered regions as GRanges
# get coordinates from FPKM table
coverage_gr <- 
  fpkm_tb %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.)

# sort
coverage_gr <- sortSeqlevels(coverage_gr)
coverage_gr <- sort(coverage_gr)


### set parameters
mirna_name <- "ssc-mir-205"
mirna_seq <- "uccuucauuccaccggagucug"

# get 4 patterns of one miRNA - 8mer, 7mer-m8 and 7mer-1a
mir_pattern <- 
  mirna_seq %>% 
  toupper(.) %>% 
  dplyr::tibble(mature_sequence = .) %>% 
  dplyr::mutate(seed_6mer = str_sub(mature_sequence, 2, 7) %>% str_replace_all(., "U", "T"), 
                seed_7mer_m8 = str_sub(mature_sequence, 2, 8) %>% str_replace_all(., "U", "T"), 
                seed_8mer = str_c("T", str_sub(mature_sequence, 2, 8) %>% str_replace_all(., "U", "T")), 
                seed_7mer_1a = str_c("T", str_sub(mature_sequence, 2, 7) %>% str_replace_all(., "U", "T"))) %$%
  DNAStringSet(c(seed_6mer, seed_7mer_m8, seed_8mer, seed_7mer_1a)) %>% 
  reverseComplement(.) %>% 
  magrittr::set_names(c("6mer", "7mer_m8", "8mer", "7mer_1a")) %>% 
  .[names(.) %in% c("6mer", "7mer_m8")]


# look for seed patterns in  
mirna_targets_list <- purrr::map(names(mir_pattern), function(mir_seed_name){
  
  # find one pattern
  seed_match_gr <- Biostrings::vmatchPattern(pattern = mir_pattern[[mir_seed_name]], 
                                             subject = BSgenome.Sscrofa.UCSC.susScr11, 
                                             max.mismatch = 0)
  
  # sort
  seed_match_gr <- sortSeqlevels(seed_match_gr)
  seed_match_gr <- sort(seed_match_gr, ignore.strand = T)
  
  # overlap each pattern with coverage regions
  overlaps <- findOverlaps(seed_match_gr, coverage_gr, ignore.strand = T)
  
  # extract seeds and their respective coverage regions, summarize
  seed_match_tb <- 
    tibble(gene_id = mcols(coverage_gr[subjectHits(overlaps)])$gene_id, 
           hit_strand = seed_match_gr[queryHits(overlaps)] %>% strand(.) %>% as.character(.)) %>% 
    dplyr::group_by(gene_id) %>% 
    dplyr::summarise(count = n(), 
                     hits = str_c(hit_strand, collapse = "")) %>% 
    dplyr::ungroup(.) %>% 
    set_names(., c("gene_id", str_c("count_", mir_seed_name), str_c("hits_", mir_seed_name)))
  
  # return
  return(seed_match_tb)
  
})

# join to one table
mirna_targets_tb <- 
  mirna_targets_list %>% 
  purrr::reduce(., full_join, by = "gene_id") %>% 
  dplyr::mutate_all(~(replace(., is.na(.), 0))) %>% 
  dplyr::arrange(desc(`+`(count_6mer, count_7mer_m8)))

# add annotation and FPKM 
mirna_targets_fpkm <- 
  mirna_targets_tb %>% 
  dplyr::left_join(., fpkm_tb, by = "gene_id") %>% 
  dplyr::select(coordinates, FPKM_in_GV = GV, count_6mer, count_7mer_m8, hits_6mer, hits_7mer_m8) 

# save
readr::write_csv(mirna_targets_fpkm, file.path(outpath, str_c("Sscrofa11.1.20200424.genome", mirna_name, "targets.csv", sep = ".")))
