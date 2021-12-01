### INFO: 
### DATE: Thu Jul 18 19:52:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/miR-205_pig")

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


### exon coordinates
# genome path
genome_dir <- "/common/DB/genome_reference/pig/susScr11.Sscrofa11.1.GCA_000003025.6"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames\\.geneInfo\\.csv$"), full.names = T)

# exon coordinates path
exon_gr_path <- list.files(genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*\\.reducedExons\\.RDS"), full.names = T)


### assembled transcriptome
# assembly path
assembly_path <- "/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/transcriptome_assemblies/pig.susScr11/trinity.genome_guided"

# transcript sequences path
trancripts_seq_path <- list.files(path = assembly_path, pattern = str_c(".*fasta$"), full.names = T)

# mapped transcripts path
transcripts_bam_path <- list.files(path = assembly_path, pattern = str_c(".*bam$"), full.names = T)


### FPKM
# maternal transcriptome path
maternal_transcriptome_path <- "/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/ensembl_counts.featureCounts/pig.susScr11"

# FPKM table path
fpkm_tb_path <- list.files(path = maternal_transcriptome_path, pattern = str_c("ensembl\\.", ensembl_version, ".*FPKM_mean\\.csv$"), full.names = T)

######################################################## READ DATA
# read transcripts info
genes_info <- readr::read_csv(genes_info_path)

# read exon coordinates
exons_gr <- readRDS(exon_gr_path)

# read transcript sequences
transcripts_seq <- Biostrings::readDNAStringSet(trancripts_seq_path)

# read transcripts mapping
transcripts_bam <- GenomicAlignments::readGAlignments(transcripts_bam_path, use.names = T)

# read FPKM table
fpkm_tb <- readr::read_csv(fpkm_tb_path)

######################################################## MAIN CODE
# overlap transcripts with exons
overlaps <- findOverlaps(transcripts_bam, unlist(exons_gr))


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
  magrittr::set_names(c("6mer", "7mer-m8", "8mer", "7mer-1a"))

# look for seed patterns in 3'UTRs 
mirna_targets_list <- purrr::map(names(mir_pattern), function(mir_seed_name){
  
  # find one pattern
  seed_match_tb <- 
    Biostrings::vmatchPattern(pattern = mir_pattern[[mir_seed_name]], 
                              subject = transcripts_seq, 
                              max.mismatch = 0) %>% 
    unlist(.) %>% 
    as_tibble(.) %>% 
    dplyr::select(-width) %>% 
    dplyr::group_by(names) %>% 
    dplyr::summarize(count = n()) %>% 
    set_names(., c("transcript_name", mir_seed_name))

})

# join to one table
mirna_targets_tb <- 
  mirna_targets_list %>% 
  purrr::reduce(., full_join, by = "transcript_name") %>% 
  dplyr::mutate_all(~(replace(., is.na(.), 0))) %>% 
  dplyr::arrange(desc(`+`(`6mer`, `7mer-m8`)))

# add annotation and FPKM 
mirna_targets_fpkm <- 
  mirna_targets_tb %>% 
  dplyr::mutate(transcript_id = str_remove(transcript_name, " .*") %>% str_remove(., "\\.[0-9]+$"), 
                gene_id = str_extract(transcript_name, "(?<=gene:)ENSSSCG.*(?= gene_biotype)") %>% str_remove(., "\\.[0-9]+$")) %>% 
  dplyr::select(gene_id, transcript_id, `6mer`:`7mer-1a`) %>% 
  dplyr::left_join(., genes_info %>% dplyr::select(gene_name, gene_id), by = "gene_id") %>% 
  dplyr::left_join(., fpkm_tb %>% dplyr::select(-gene_name), by = "gene_id") %>% 
  dplyr::select(transcript_id, gene_id, gene_name, `6mer`:`7mer-1a`, FPKM_in_GV = GV, coordinates, gene_biotype, gene_description)

# save
readr::write_csv(mirna_targets_fpkm, file.path(outpath, str_c("ensembl", ensembl_version, "Sscrofa11.1.20200424.transcriptome", mirna_name, "targets.csv", sep = ".")))
