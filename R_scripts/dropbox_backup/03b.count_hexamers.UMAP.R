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

library(umap)
library(ggrepel)

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
utrs_seq <- utrs_seq[nchar(utrs_seq) >= 200]

# divide into equally sized bins (or equally spaced bins)
fpkm_tb_bins <- 
  fpkm_tb %>% 
  dplyr::select(gene_id, fpkm = GV) %>% 
  dplyr::filter(fpkm > 0) %>%
  as.data.table(.) %>%
  .[order(fpkm)] %>%
  .[, position := order(fpkm)] %>%
  .[order(position), bin := as.numeric(Hmisc::cut2(1:.N, g = 20))] %>% # equally sized bins
  # .[order(fpkm), bin := as.numeric(Hmisc::cut2(fpkm, g = 10))] %>% # equally spaced bins
  as_tibble(.)

# count hexamers
kmer_freq <- Biostrings::oligonucleotideFrequency(utrs_seq, 7, as.prob = T, with.labels = TRUE)
rownames(kmer_freq) <- names(utrs_seq)

# get 1000 the most variable kmers
kmer_sd <- apply(kmer_freq, 2, sd)
kmer_topN <- kmer_freq[, order(kmer_sd, decreasing = T)[1:500]]
kmer_umap <- umap(kmer_topN)

# get table for plot
kmer_umap_tb <- 
  kmer_umap$layout %>% 
  as.data.frame(.) %>% 
  as_tibble(., rownames = "gene_id") %>% 
  dplyr::mutate(gene_id = str_remove(gene_id, "\\..*")) %>% 
  dplyr::filter(gene_id %in% fpkm_tb_bins$gene_id) %>% 
  dplyr::left_join(., fpkm_tb_bins, by = "gene_id") %>% 
  dplyr::mutate(bin = factor(bin, levels = 1:20))
  
# annotate
kmer_umap_labels_tb <- 
  kmer_umap_tb %>% 
  dplyr::filter(V1 < -50)

# plot
kmer_umap_plot <- 
  ggplot() + 
  geom_point(data = kmer_umap_tb, aes(V1, V2, color = bin), 
             size = 5) + 
  # geom_label_repel(data = kmer_umap_labels_tb, aes(V1, V2, label = gene_id), 
  #                  fontface = "bold", color = "black", box.padding = 0.35, 
  #                  point.padding = 0.5, segment.color = "grey50") +
  # guides(color = guide_legend(override.aes = list(shape = 23, size = 5))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = - 0.2),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom") +
  ggtitle("UMAP on 7mers") +
  ggsave(filename = file.path(outpath, "UMAP.top100_variable_hexamers.pig_FPKM_50_bins.png"), 
         width = 10, height = 12)

