### INFO: 
### DATE: Fri Jun 28 22:08:17 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Fugaku_intron_CpG/count_intronic_reads")

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

# intronic .bed path
introns_bed_path <- file.path(inpath, "../intronic_reads", "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.intronic.bed")

# summarizedOverlaps path
se_path <- file.path(inpath, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.intronic.se.RDS")

# sample table path
sample_table_path <- file.path(inpath, "../mapped_bams/log.Fugaku.stats_and_tracks.csv")

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read intronic regions from bed
introns_gr <- rtracklayer::import.bed(introns_bed_path)

# read summarizedExperiment
se <- readRDS(se_path)

# read sample table
sample_table <- readr::read_csv(sample_table_path)

######################################################## MAIN CODE
# clean sample table
sample_table %<>% 
  dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA) %>% 
  as.data.table(.)

# get unique names of introns
mcols(introns_gr)$score <- NULL
mcols(introns_gr)$name <- make.unique(mcols(introns_gr)$name)
names(introns_gr) <- mcols(introns_gr)$name
mcols(introns_gr)$name <- NULL

# get length of introns
introns_tb <- 
  introns_gr %>% 
  as.data.frame(.) %>% 
  as.data.table(., keep.rownames = "gene_id")

# FPKM
fpkm_df <-
  assay(se) %>%
  as.data.table(., keep.rownames = "gene_id") %>% 
  melt(., 
       id.vars = c("gene_id"),
       variable.name = "sample_id", 
       value.name = "counts") %>% 
  .[, sample_id := str_remove(sample_id, ".intronic.bam")] %>% 
  .[sample_table, on = "sample_id", `:=`(fpm = (counts / round(library_size / 1E6, 6)))] %>% 
  .[introns_tb, on = "gene_id", `:=`(fpkm = fpm / round(width / 1E3, 3))] %>%
  dcast(., gene_id ~ sample_id, value.var = "fpkm") %>% 
  .[introns_tb, on = "gene_id", `:=`(coordinates = str_c(i.seqnames, i.start, i.end, sep = " "), 
                                     strand = i.strand,
                                     width = i.width)] %>% 
  .[] %T>% 
  write_csv(., file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".FPKM.csv")))