#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: get expression of developmental profile
### DATE: /common/WORK/fhorvat/Projekti/Svoboda/Analyses/developmental_profile_expression/summarizedExperiments
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/developmental_profile_expression/summarizedExperiments")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gtf path
gtf_path <- list.files(path = genome_dir, pattern = "ensembl.93.*UCSCseqnames.gtf.gz$", full.names = T)

# genes info path
gene_info_path <- list.files(path = genome_dir, pattern = "ensembl.93.*UCSCseqnames.geneInfo.csv$", full.names = T)

######################################################## READ DATA
# read gtf
ensembl_gtf <- read_delim(file = gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# read gene info
gene_info_tb <- readr::read_csv(gene_info_path)

######################################################## MAIN CODE
# convert GTF to GRanges, get only exons, reduce (get only 1st and last exon for Sirena1 (ENSMUSG00000110001))
exons_gr <- 
  gtfToGRanges(ensembl_gtf, filter = "exon") %>% 
  .[!((mcols(.)$gene_id == "ENSMUSG00000110001") & !(mcols(.)$exon_id %in% c("ENSMUSE00001389346", "ENSMUSE00001389258")))] %>% 
  .[mcols(.)$gene_id == "ENSMUSG00000110001"] %>% 
  GenomicRanges::split(., .$gene_id) %>%
  GenomicRanges::reduce(., ignore.strand = T) %>%
  unlist(.)

# add strand info, split again
exons_gr$gene_id <- names(exons_gr)
names(exons_gr) <- NULL
exons_gr <-
  exons_gr %>%
  as.data.frame(.) %>%
  as_tibble(.) %>%
  dplyr::select(-strand) %>%
  dplyr::left_join(., gene_info_tb %>% dplyr::select(gene_id, strand), by = "gene_id") %>%
  GenomicRanges::GRanges(.) %>%
  GenomicRanges::split(., .$gene_id)

# save RDS
saveRDS(object = exons_gr, file = file.path(outpath, gtf_path %>% basename(.) %>% stringr::str_replace(., ".gtf.gz", ".only_short_Sirena1.reducedExons.RDS")))
