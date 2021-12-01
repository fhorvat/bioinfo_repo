### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/siRNA.Valeria/datasets/mouse_mESC.Dicer_mutants.small_RNAseq.2021_May/Analysis/Dicer_KO_reads")

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
library(GenomicAlignments)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames\\.geneInfo\\.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames\\.reducedExons\\.RDS$"), full.names = T)

# rmsk path
rmsk_path <- list.files(path = genome_dir, pattern = "rmsk.*\\.clean\\.fa\\.out\\.gz", full.names = T)

# miRBase path
mirbase_path <- list.files(path = genome_dir, pattern = "miRBase.*\\gff3", full.names = T)

# list bam files
bam_list <- list.files(path = inpath, pattern = ".*\\.21to23nt\\.bam", full.names = T)

######################################################## READ DATA

######################################################## MAIN CODE
# read gene info
genes_info <- readr::read_csv(genes_info_path)

# read exons
exons_gr <-
  readRDS(file = exons_path) %>%
  tibble::as_tibble(.) %>%
  dplyr::left_join(., genes_info %>% dplyr::select(gene_id, gene_biotype), by = "gene_id") %>%
  dplyr::select(seqnames:strand, gene_id, gene_biotype) %>%
  GenomicRanges::GRanges(.)

# remove miRNA
exons_gr <- exons_gr[mcols(exons_gr)$gene_biotype != "miRNA"]

# read repeatMasker
rmsk_gr <-
  read_delim(file = rmsk_path, delim = "\t", col_types = cols(start = col_double(), end = col_double())) %>%
  dplyr::mutate(gene_id = str_c(seqnames, ":", start, "-", end, "|", strand, "|", repName)) %>%
  dplyr::select(seqnames:strand, gene_id, gene_biotype = repClass) %>%
  GenomicRanges::GRanges(.)

# read and clean miRBase miRNA annoation
mirna_gr <- rtracklayer::import.gff(con = mirbase_path)
mcols(mirna_gr) <- mcols(mirna_gr)[, c("ID")]
names(mcols(mirna_gr)) <- "gene_id"
mcols(mirna_gr)$gene_biotype <- "miRNA.mature"

# join repeatMasker, exons and miRNA annotation
features_gr <- c(mirna_gr, rmsk_gr, exons_gr)


### loop through bam files
purrr::map(bam_list, function(bamfile)){
  
  # read bam
  chunk <- readGAlignmentsList(file = bamfile,
                               param = ScanBamParam(what = "qname"))
  
  # unlist, set names of reads
  chunk <- unlist(chunk)
  names(chunk) <- mcols(chunk)$qname
  
  # transform to grglist (which gets ranges of only alignment part)
  chunk <-
    GenomicRanges::grglist(chunk) %>%
    unlist(.)
  
  # find overlaps between two GRanges - chunk of alignments and annotation
  hits <- findOverlaps(chunk, features_gr, ignore.strand = T, type = "any")

  # get hits in alignments
  read_hits <-
    extractList(chunk, as(queryHits(hits), "List")) %>%
    unlist(.)

  # get hits in annotation
  subject_hits <-
    extractList(features_gr, as(subjectHits(hits), "List")) %>%
    unlist(.)

  # create a table
  hits_tb <-
    tibble(read_names = names(read_hits),
           hit_biotype = mcols(subject_hits)$gene_biotype) %>%
    dplyr::group_by(read_names) %>%
    dplyr::summarise(n_biotype = length(unique(hit_biotype)),
                     hit_biotype = str_c(hit_biotype, collapse = ", "))
  
  counts <- countOverlaps(features_gr, chunk, ignore.strand = F, type = "any")
  
  features_gr$counts <- counts
  
  counts_tb <- 
    features_gr %>% 
    as_tibble(.) %>% 
    dplyr::filter(counts > 0) %>% 
    dplyr::arrange(-counts)
  
}






