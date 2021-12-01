### INFO: 
### DATE: Tue Oct 29 22:03:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/tmp/test")

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
library(GenomicFeatures)
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
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)

# .gtf path
gtf_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.gtf.gz$"), full.names = T)

# RDS path
rds_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis/ensembl.GRCm38.89.CNOT6L.summarizedOverlaps.RDS"

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read .gtf
gtf_gr <- rtracklayer::import.gff(gtf_path)

# read assay from .RDS
se <- readRDS(rds_path)

######################################################## MAIN CODE
# get lengths of all exons for each gene on chromosome 5
exon_lengths <-
  gtf_gr[seqnames(gtf_gr) %in% c("chr5")] %>% 
  as_tibble(.) %>% 
  dplyr::filter(type == "exon") %>%
  dplyr::select(seqnames, start, end, width, strand, gene_id, transcript_id, exon_id)

# filter se
se_filt <- se[rownames(se) %in% exon_lengths$gene_id, ]
# colnames(se_filt) <- str_remove(colnames(se_filt, "\\.Aligned\\.sortedByCoord\\.out\\.bam"))
# saveRDS(se_filt, file = file.path(outpath, "ensembl.93.GRCm38.p6.CNOT6L.summarizedOverlaps.chr5.RDS"))

# save
assay(se_filt) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  as_tibble(.) %>%
  dplyr::rename_with(.fn = ~str_remove(string = .x, pattern = "\\.Aligned\\.sortedByCoord\\.out\\.bam"),
                     .cols = starts_with("s_")) %>%
  write_csv(., "ensembl.93.GRCm38.p6.CNOT6L.summarizedOverlaps.chr5.csv")


### version 1 - save exon lengths as it is, they need to find longest splice variant for each gene
# save 
readr::write_csv(exon_lengths, file = file.path(outpath, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.exon_lengths.chr5.csv"))


### version 2 - make the exon lengths table wide so they have to tidy it first, but filter longest splice variant 
# get the longest transcript for each gene
transcript_lengths <-
  exon_lengths %>%
  dplyr::group_by(gene_id, transcript_id) %>%
  dplyr::summarise(width = sum(width)) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::slice_max(order_by = width, n = 1, with_ties = F) %>%
  dplyr::ungroup(.)

# filter exon lengths
exon_lengths %<>%
  dplyr::filter(transcript_id %in% transcript_lengths$transcript_id) %>%
  dplyr::select(seqnames, start, end, strand, gene_id, exon_id)

# pivot to wide table
exon_lengths_wide <- 
  exon_lengths %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::mutate(exon_index = ifelse(strand =="+", 1:n(), n():1), 
                exon_index = str_c("exon", exon_index)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::pivot_wider(id_cols = c(seqnames, strand, gene_id), 
                     names_from = exon_index, 
                     values_from = c(start, end), 
                     names_sep  = "_")

# # how to get tidy table
# exon_lengths_wide %>%
#   tidyr::pivot_longer(cols = starts_with(c("start", "end")),
#                       names_to = c("start_or_end", "exon_id"),
#                       names_sep = "_",
#                       values_to = c("coordinate"),
#                       values_drop_na = T) %>%
#   tidyr::pivot_wider(id_cols = c(seqnames:gene_id, exon_id),
#                      values_from = coordinate,
#                      names_from = start_or_end)

# save exon lengths 
readr::write_csv(exon_lengths_wide, file = file.path(outpath, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.exon_lengths_wide.chr5.csv"))


### version 3 - filter gtf to include only genes from 5th chromosome, require to read in .gtf without rtracklayer 
# and find the longest transcript for each gene
ensembl_gtf <- read_delim(file = gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))
gtf_chr5 <- 
  ensembl_gtf %>% 
  dplyr::filter(X1 == "chr5")

write.table(x = gtf_chr5, file = file.path(outpath, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.chr5.gtf"), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

