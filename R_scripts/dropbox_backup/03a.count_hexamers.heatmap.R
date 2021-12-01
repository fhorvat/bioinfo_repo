### INFO: 
### DATE: Sat Apr 18 22:14:36 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/2020_paper/miRNA_targets_expression")

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

library(biomaRt)
library(Biostrings)

library(pheatmap)
library(viridis)
library(dendsort)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# sort clusters
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

######################################################## PATH VARIABLES
# set inpath 
inpath <- file.path(getwd(), "sequences")

# set outpath
outpath <- file.path(getwd(), "sequences")

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/pig/susScr11.Sscrofa11.1.GCA_000003025.6"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# get 3' UTR fasta path
utrs_path <- file.path(inpath, "ensembl.93.Sscrofa11.1.20200423.3UTR.fasta")

# FPKM expression path
fpkm_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/ensembl_counts.featureCounts/pig.susScr11", 
                       "ensembl.93.Sscrofa11.1.20180920.UCSCseqnames.FPKM_mean.csv")

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read 3' UTRs
utrs_seq <- Biostrings::readDNAStringSet(utrs_path)

# read FPKM table
fpkm_tb <- readr::read_csv(fpkm_path)

######################################################## MAIN CODE
# remove UTRs shorter than 20 nt
utrs_seq <- utrs_seq[nchar(utrs_seq) >= 20]

# count hexamers
kmer_freq <- Biostrings::oligonucleotideFrequency(utrs_seq, 7, as.prob = F, with.labels = TRUE)

# divide into equally sized bins (or equally spaced bins)
fpkm_tb_bins <- 
  fpkm_tb %>% 
  dplyr::select(gene_id, fpkm = GV) %>% 
  # dplyr::filter(fpkm > 0) %>%
  as.data.table(.) %>%
  .[order(fpkm)] %>%
  .[, position := order(fpkm)] %>%
  # .[order(position), bin := as.numeric(Hmisc::cut2(1:.N, g = 50))] %>% # equally sized bins
  .[order(fpkm), bin := as.numeric(Hmisc::cut2(fpkm, g = 50))] %>% # equally spaced bins
  as_tibble(.)

# sum hexamer frequency by bins
kmer_freqs_tb <- 
  kmer_freq %>% 
  as_tibble(.) %>% 
  dplyr::mutate(gene_id = names(utrs_seq) %>% str_remove(., "\\..*")) %>% 
  dplyr::select(gene_id, everything()) %>% 
  dplyr::filter(gene_id %in% fpkm_tb_bins$gene_id) %>% 
  dplyr::left_join(., fpkm_tb_bins, by = "gene_id") %>%
  dplyr::select(-c(gene_id, position, fpkm)) %>%
  tidyr::pivot_longer(., cols = -bin, names_to = "kmer", values_to = "count") %>% 
  dplyr::filter(!(kmer %in% c("AAAAAAA", "TTTTTTT"))) %>%
  # dplyr::filter(!str_detect(kmer, "[A,T]{6}")) %>%
  dplyr::group_by(bin) %>% 
  dplyr::mutate(kmer_sum = sum(count)) %>% 
  dplyr::group_by(bin, kmer) %>% 
  dplyr::summarise(count = sum(count), kmer_sum = unique(kmer_sum)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(kmer_freq = count / kmer_sum) %>% 
  dplyr::select(bin, kmer, kmer_freq) %>% 
  tidyr::pivot_wider(., id_cols = bin, names_from = kmer, values_from = kmer_freq)

# to matrix for heatmap
kmer_freqs_mat <- 
  kmer_freqs_tb %>% 
  as.data.frame(.) %>%
  column_to_rownames(., var = "bin") %>%
  as.matrix(.) %>%
  t(.)

# filter top N rows by variance
ntop <- 50
rv <- genefilter::rowVars(kmer_freqs_mat) 
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
kmer_freqs_mat_filt <- kmer_freqs_mat[select, ]

# sort columns
mat_cluster_rows <- sort_hclust(hclust(dist(kmer_freqs_mat_filt)))

# plot
pheatmap::pheatmap(kmer_freqs_mat_filt,
                   col = viridis(50),
                   cluster_cols = F,
                   cluster_rows = T,
                   show_rownames = T, 
                   show_colnames = T, 
                   border_color = NA, 
                   file = file.path(outpath, "heatmap.top50_variable_heptamers.pig_FPKM_50_bins.png"),
                   height = 20,
                   width = 20)
