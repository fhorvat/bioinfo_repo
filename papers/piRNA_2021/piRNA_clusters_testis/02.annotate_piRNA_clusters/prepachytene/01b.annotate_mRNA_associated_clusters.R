### INFO: 
### DATE: Mon Aug 31 14:57:45 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/annotate_clusters")

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

library(openxlsx)
library(rtracklayer)
library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# clusters path
clusters_path <- file.path(inpath, "MesAur1.1k_windows.rpkm_cutoff.1_supplementary_pre-pachytene.xlsx")

# genome path
genome_path <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"

# gtf path
gtf_path <- file.path(genome_path, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.gtf")

# gene info path
gene_info_path <- str_replace(gtf_path, "\\.gtf$", ".geneInfo.csv")
  
# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq/Data/Mapped/STAR_mesAur1.new/8_merged_replicates"

# WT coverage path
wt_coverage_path <- file.path(mapped_path, "s_testis_Mov10l_WT_13dpp.24to31nt.scaled.bw")

# testis RNA-seq expression FPKM table path
testis_fpkm_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.RNAseq/Analysis/expression"
testis_fpkm_path <- list.files(testis_fpkm_path, ".*\\.FPKM_mean\\.csv$", full.names = T)

# testis small RNA-seq expression path
testis_fpm_small_path_sense <- file.path(inpath, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.sense.FPM.csv")
testis_fpm_small_path_antisense <- file.path(inpath, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.antisense.FPM.csv")

######################################################## READ DATA
# read clusters table
clusters_tb <- 
  openxlsx::read.xlsx(clusters_path) %>% 
  as_tibble(.)

# read gtf
gtf_gr <- rtracklayer::import.gff(gtf_path)

# read gene info
gene_info <- readr::read_csv(gene_info_path)

# read WT coverage
wt_coverage <- rtracklayer::import.bw(wt_coverage_path)

# read RNA-seq expression in testis
testis_fpkm <- readr::read_csv(testis_fpkm_path)

# read FPM tables
testis_fpm_small_sense <- readr::read_csv(testis_fpm_small_path_sense)
testis_fpm_small_antisense <- readr::read_csv(testis_fpm_small_path_antisense)

######################################################## MAIN CODE
### tidy and prepare data
# clean clusters table, get GRanges
clusters_gr <- 
  clusters_tb %>% 
  dplyr::select(coordinates) %>% 
  tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ", remove = F) %>% 
  GenomicRanges::GRanges(.)

# get gene coordinates
gene_gr <- gtf_gr[mcols(gtf_gr)$type == "gene"]

# get exon coordinates
exon_gr <- gtf_gr[mcols(gtf_gr)$type == "exon"]

# get 3'UTR coordinates
utr_gr <- gtf_gr[mcols(gtf_gr)$type == "three_prime_utr"]

# get non-zero coverage
wt_coverage_expressed <- 
  wt_coverage[wt_coverage$score > 0] %>% 
  reduce(.)

# clean sense/antisense small RNA-seq FPM
testis_fpm_small_sense %<>% 
  dplyr::select(gene_id, 
                small_RPM.sense.Mov10l_WT_13dpp = s_testis_Mov10l_WT_13dpp, 
                small_RPM.sense.Mov10l_KO_13dpp = s_testis_Mov10l_KO_13dpp, 
                small_RPM.sense.Mov10l_WT_21dpp = s_testis_Mov10l_WT_21dpp,
                small_RPM.sense.Mov10l_KO_21dpp = s_testis_Mov10l_KO_21dpp)

# clean sense/antisense small RNA-seq FPM
testis_fpm_small_antisense %<>% 
  dplyr::select(gene_id, 
                small_RPM.antisense.Mov10l_WT_13dpp = s_testis_Mov10l_WT_13dpp, 
                small_RPM.antisense.Mov10l_KO_13dpp = s_testis_Mov10l_KO_13dpp, 
                small_RPM.antisense.Mov10l_WT_21dpp = s_testis_Mov10l_WT_21dpp,
                small_RPM.antisense.Mov10l_KO_21dpp = s_testis_Mov10l_KO_21dpp)


### get clusters overlapping with genes
# overlap
overlaps <- GenomicRanges::findOverlaps(clusters_gr, gene_gr)

# get cluster/genes overlaps
clusters_genes_tb <- tibble(coordinates = clusters_gr[queryHits(overlaps)]$coordinates,
                            gene_id = gene_gr[subjectHits(overlaps)]$gene_id)

# get exons which have coverage in 13 dpp WT small RNA-seq
exon_gr_filt <- 
  exon_gr[exon_gr$gene_id %in% clusters_genes_tb$gene_id] %>% 
  subsetByOverlaps(., wt_coverage_expressed) %$% 
  gene_id %>% 
  unique(.) %>% 
  tibble(gene_id = ., exon_match = "Y")

# get 3'UTRs which have coverage in 13 dpp WT small RNA-seq
utr_gr_filt <- 
  utr_gr[utr_gr$gene_id %in% clusters_genes_tb$gene_id] %>% 
  subsetByOverlaps(., wt_coverage_expressed) %$% 
  gene_id %>% 
  unique(.) %>% 
  tibble(gene_id = ., `3pUTR_match` = "Y") 

# add info about exon/3pUTR coverage to gene
clusters_annotated <- 
  clusters_genes_tb %>% 
  dplyr::left_join(., exon_gr_filt, by = "gene_id") %>% 
  dplyr::left_join(., utr_gr_filt, by = "gene_id") %>% 
  dplyr::mutate(exon_match = replace(exon_match, is.na(exon_match), "N"), 
                `3pUTR_match` = replace(`3pUTR_match`, is.na(`3pUTR_match`), "N")) %>% 
  dplyr::left_join(., gene_info %>% dplyr::select(gene_id, gene_name), by = "gene_id") %>% 
  dplyr::group_by(coordinates) %>% 
  dplyr::summarise(gene_name = str_c(unique(gene_name), collapse = ","), 
                   gene_id = str_c(unique(gene_id), collapse = ","), 
                   exon_match = any(exon_match == "Y"), 
                   `3pUTR_match` = any(`3pUTR_match` == "Y")) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(exon_match | `3pUTR_match`) %>% 
  dplyr::select(coordinates, exon_match, `3pUTR_match`, gene_name, gene_id)

# add to the original cluster table
clusters_tb_final <- 
  clusters_tb %>% 
  dplyr::left_join(., clusters_annotated, by = "coordinates") %>% 
  dplyr::mutate(exon_match = replace(exon_match, is.na(exon_match), F), 
                `3pUTR_match` = replace(`3pUTR_match`, is.na(`3pUTR_match`), F))

# # save
# openxlsx::write.xlsx(clusters_tb_final, file = file.path(outpath, "MesAur1.1k_windows.rpkm_cutoff.1_supplementary_pre-pachytene.genes_annotated.20200831.FH.xlsx"))


### get tables with piRNA associated genes
# get all genes in table
gene_pirna_tb <- 
  clusters_tb_final %$% 
  gene_id %>% 
  str_split(., ",") %>% 
  unlist(.) %>% 
  .[!is.na(.)] %>% 
  tibble(gene_id = .)

# add info to the table
gene_pirna_info <- 
  gene_pirna_tb %>% 
  dplyr::left_join(., gene_info, by = "gene_id") %>% 
  tidyr::unite(col = "coordinates", seqnames, start, end, sep = " ") %>% 
  dplyr::select(gene_id, gene_name, coordinates, strand) %>% 
  dplyr::left_join(., testis_fpkm %>% dplyr::select(gene_id, FPKM.Mov10l_WT_13dpp = Mov10l_WT_13dpp, FPKM.Mov10l_WT_21dpp = Mov10l_WT_21dpp), by = "gene_id") %>% 
  dplyr::left_join(., testis_fpm_small_sense, by = "gene_id") %>% 
  dplyr::left_join(., testis_fpm_small_antisense, by = "gene_id") %>% 
  dplyr::select(gene_id:FPKM.Mov10l_WT_21dpp,
                small_RPM.sense.Mov10l_WT_13dpp, small_RPM.antisense.Mov10l_WT_13dpp, 
                small_RPM.sense.Mov10l_KO_13dpp, small_RPM.antisense.Mov10l_KO_13dpp,
                small_RPM.sense.Mov10l_WT_21dpp, small_RPM.antisense.Mov10l_WT_21dpp, 
                small_RPM.sense.Mov10l_KO_21dpp, small_RPM.antisense.Mov10l_KO_21dpp)

# write
openxlsx::write.xlsx(gene_pirna_info, file = file.path(outpath, "piRNA_associated_mRNAs.MesAur1.1k_windows.rpkm_cutoff.1_supplementary_pre-pachytene.20200831.FH.xlsx"))

  
