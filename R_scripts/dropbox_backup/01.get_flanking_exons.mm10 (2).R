### INFO: 
### DATE: Wed Jul 10 00:18:10 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/lnc1_locus/other/flanking_genes.mRNA")

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
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)
library(seqinr)

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

# .gtf path
gtf_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.gtf.gz$"), full.names = T)

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read .gtf
gtf_tb <- read_delim(file = gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

######################################################## MAIN CODE
# convert GTF to GRanges, get only exons
exons_gr <-
  gtfToGRanges(gtf_tb, filter = "exon") %>%
  GenomicRanges::split(., .$gene_id)

### get genes of interests:
# Uevld  = ENSMUSG00000043262, downstream, exon 1, minus strand (same strand as lnc1)
# Spty2d1 = ENSMUSG00000049516, upstream, exon 6, minus strand (same strand as lnc1)
genes_list <- c("ENSMUSG00000043262", "ENSMUSG00000049516")

# filter exons list
exons_filt <- 
  exons_gr[names(exons_gr) %in% genes_list] %>% 
  unlist(.) %>% 
  as_tibble(.) %>% 
  group_by(gene_id) %>% 
  arrange(start) %>% 
  dplyr::mutate(exon_total = n(),
                exon_n = n():1) %>% 
  dplyr::ungroup(.) 

# for ENSMUSG00000043262 take exon 2 (first exon conserved in rat)
Uevld_exon <- 
  exons_filt %>% 
  dplyr::filter(gene_id == "ENSMUSG00000043262", exon_n == 2) %>% 
  GRanges(.) %T>%
  export.bed(., file.path(outpath, "Uevld.exon2.bed"))

# for ENSMUSG00000049516 take last exon
Spty2d1_exon <- 
  exons_filt %>% 
  dplyr::filter(gene_id == "ENSMUSG00000049516", exon_n == exon_total) %>% 
  GRanges(.) %T>%
  export.bed(., file.path(outpath, "Spty2d1.exon7.bed"))

# get sequences

