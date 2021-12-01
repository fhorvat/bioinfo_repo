### INFO: 
### DATE: Mon Oct 29 16:48:51 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/miRNA_binding_sites/last_exon.let7_counts")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(ggthemes)

library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "gtfToGRanges.R"))

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# ensembl .gtf path
ensembl_path <- list.files(path = genome_path, pattern = "ensembl.93.*[0-9]{6}.UCSC.*gtf.gz", full.names = T)

# gene info path
gene_info_path <- list.files(path = genome_path, pattern = "ensembl.93.*[0-9]{6}.UCSC.*geneInfo.csv", full.names = T)

# oocyte expressed genes
oocyte_fpkm_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/miRNA_binding_sites/last_exon.let7_counts/ensembl.93.mouse_B6.mm10.avgFPKM_above_0.1.RDS"

######################################################## READ DATA
# read ensembl .gtf
ensembl_gtf <- read_delim(file = ensembl_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# read gene info
genes_info <- readr::read_csv(gene_info_path)

# read oocyte expressed genes
oocyte_fpkm <- readRDS(oocyte_fpkm_path)

######################################################## MAIN CODE
# get all exons from gtf
ensembl_last_exons_gr <-
  ensembl_gtf %>%
  gtfToGRanges(., filter = "exon")

# get last exon for each gene
ensembl_last_exons <- 
  ensembl_last_exons_gr %>% 
  as.data.table(.) %>% 
  .[gene_id %in% oocyte_fpkm$gene_id] %>% 
  .[order(start), .SD, .(gene_id)] %>%
  .[strand == "+" , `:=`(exon_count = .N , exon_idx = 1:.N) , by = gene_id[strand == "+"]] %>% 
  .[strand == "-" , `:=`(exon_count = .N , exon_idx = .N:1) , by = gene_id[strand == "-"]] %>%
  .[exon_idx == exon_count] %>%
  .[] %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T)

# get sequences of last exon in the gene
ensembl_last_exons_seq <- BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, ensembl_last_exons)
names(ensembl_last_exons_seq) <- ensembl_last_exons$gene_id

### let-7 targets
# miRNA name
mirna_name <- "let7"

# miRNA sequence 
mirna_seq <- 
  "UGAGGUAGUAGGUUGUAUAGUU" %>% 
  str_replace_all(., "U", "T") %>% 
  Biostrings::DNAStringSet(.) %>% 
  unlist(.)

# get pattern
mirna_seed <- 
  mirna_seq %>% 
  subseq(., start = 2, end = 7) %>% 
  reverseComplement(.)

# find last exon matching pattern
match_2to7 <- 
  Biostrings::vmatchPattern(pattern = mirna_seed, 
                            subject = ensembl_last_exons_seq, 
                            max.mismatch = 0) %>% 
  unlist(.) %>% 
  as.tibble(.) %>% 
  dplyr::rename(gene_id = names) %>%
  dplyr::count(gene_id, sort = T) %>% 
  dplyr::left_join(., oocyte_fpkm, by = "gene_id") %>% 
  dplyr::left_join(., genes_info, by = "gene_id") %>% 
  dplyr::rename(seed_match_count = n, 
                CNOT6L.GV_WT.FPKM = GV_WT) %T>% 
  write_csv(., path = file.path(outpath, str_c(mirna_name, 
                                               "2to7_seed",
                                               mirna_seed %>% reverseComplement(.) %>% as.character(.),
                                               "3pExons",
                                               "CNOT6L_GV_expressed.0.1_FPKM_cut",
                                               "20190425.csv", 
                                               sep = ".")))
