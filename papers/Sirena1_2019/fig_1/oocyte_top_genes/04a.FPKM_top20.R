#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: 
### DATE: Thu Apr 25 16:16:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY

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

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### experiment
# set experiment name
experiment <- "Fugaku"

# set grouping variable(s)
grouping_variables <- "stage"


### working dir
# set working directory 
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/Fugaku_FPKM")


### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### other experiment paths
# set base experiment path
base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq", experiment)

# mapped path
mapped_path <- file.path(base_path, "Data/Mapped/STAR_mm10")

# documentation path
documentation_path <- file.path(base_path, "Data/Documentation")

# analysis path
analysis_path <- file.path(base_path, "Analysis/expression")


### documentation
# set ensembl version
ensembl_version <- 93

# sample table path
sample_table_path <- list.files(path = documentation_path, pattern = ".*sampleTable\\.csv", full.names = T)

# stats and tracks path
stats_and_tracks_path <- list.files(path = mapped_path, pattern = ".*\\.stats_and_tracks\\.csv", full.names = T)

# summarizedExperiment path
se_path <- list.files(path = analysis_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.*.se.RDS$"), full.names = T)

# FPKM table path
fpkm_path <- list.files(path = analysis_path, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames.*\\.FPKM_long\\.csv$"), full.names = T)


### genome
# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)

# repeatMasker path
rmsk_path <- list.files(path = genome_dir, pattern = "rmsk\\..*clean.*.gz$", full.names = T)


######################################################## READ DATA
# FPKM long table
fpkm_long_tb <- readr::read_csv(fpkm_path)

# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read repeatMasker
rmsk_tb <- read_delim(file = rmsk_path, delim = "\t", col_types = cols(start = col_double(), end = col_double()))

######################################################## MAIN CODE
### get gene IDs of genes whose exons overlaping rRNA
# get rRNA genes from repeatMasker, transform to GRanges
rmsk_rRNA <- 
  rmsk_tb %>% 
  dplyr::filter(repClass == "rRNA") %>% 
  GRanges(.)

# get genes which overlap annotated rRNA
genes_rRNA <- 
  exons_gr %>% 
  findOverlaps(., rmsk_rRNA, ignore.strand = T) %>% 
  queryHits(.) %>% 
  exons_gr[.] %>% 
  names(.)


### get lnc1 width, only exons 1 and 4
lnc1_exons <- 
  exons_gr[names(exons_gr) == "ENSMUSG00000110001"] %>% 
  as_tibble(.) %>% 
  dplyr::arrange(start) %>% 
  dplyr::mutate(exon_n = n():1) %>% 
  dplyr::group_by(gene_id) %>% 
  summarize(width_1_4 = sum(width[exon_n %in% c(1, 4)]))


### get top 20 genes by expression in GV.WE from Fugaku's dataset
# get top 20 genes by FPKM
fpkm_top20 <- 
  fpkm_long_tb %>%
  dplyr::left_join(., genes_info, by = "gene_id") %>% 
  dplyr::filter(sample_id == "s_GV.WE.PE", 
                !(gene_id %in% genes_rRNA), 
                (gene_biotype == "protein_coding") | (gene_id == "ENSMUSG00000110001")) %>% 
  arrange(desc(fpkm)) %>% 
  dplyr::top_n(20, fpkm)

# calculate FPKM for lnc1 using only exons 1 and 4
fpkm_lnc1 <- 
  fpkm_top20 %>% 
  dplyr::filter(gene_id == "ENSMUSG00000110001") %>% 
  dplyr::left_join(., lnc1_exons, by = "gene_id") %>% 
  dplyr::mutate(fpkm = (fpm / round(width_1_4 / 1E3, 3))) %>% 
  dplyr::mutate(gene_id = str_c(gene_id, ".exons_1_4"), 
                gene_name = str_c(gene_name, ".exons_1_4")) %>% 
  dplyr::select(-width_1_4)

# join with long table, sort, save
fpkm_final <- 
  rbind(fpkm_top20, fpkm_lnc1) %>% 
  tidyr::unite(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::select(gene_id, gene_name,	fpkm, coordinates, strand, gene_biotype, gene_description, sample_id) %>% 
  arrange(desc(fpkm))

# save long FPKM table
fpkm_final %>% 
  readr::write_csv(., file.path(outpath, basename(fpkm_path) %>% str_replace(., ".FPKM_long.csv", ".GV.top20.FPKM.csv")))
