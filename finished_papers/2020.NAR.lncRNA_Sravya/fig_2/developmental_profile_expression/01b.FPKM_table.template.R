#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: 
### DATE: Thu Apr 25 16:16:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/developmental_profile_expression/FPKMs")

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
standard.error <- function(x) {
  sqrt(var(x) / length(x))
}

######################################################## PATH VARIABLES
### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### experiment
# set ensembl version
ensembl_version <- 93

# set experiment name
experiment <- "%EXPERIMENT"

# set base experiment path
base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq", experiment)

# documentation path
documentation_path <- file.path(base_path, "Data/Documentation")

# mapped path
mapped_path <- file.path(base_path, "Data/Mapped/STAR_mm10")

# analysis path
se_main_path <- file.path(inpath, "..", "summarizedExperiments")

# summarizedExperiment path
se_path <- list.files(path = se_main_path, pattern = str_c("ensembl.", ensembl_version, ".*", "\\.short_Sirena1\\.", experiment, "\\.se\\.RDS"), full.names = T)


### documentation
# sample table path
sample_table_path <- list.files(documentation_path, ".*sampleTable.csv", full.names = T)

# reads stats path
reads_stats_path <- list.files(mapped_path, ".*\\.stats_and_tracks\\.csv", full.names = T)


### genome
# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- file.path(se_main_path, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.short_Sirena1.reducedExons.RDS")

######################################################## READ DATA
# read sample table
sample_table <- data.table::fread(sample_table_path)

# read stats and tracks table
reads_stats <- data.table::fread(reads_stats_path)

# summarizedExperiment 
se <- readRDS(se_path)

# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

######################################################## MAIN CODE
# join sample table with stats and tracks
sample_table[reads_stats, on = "sample_id", `:=`(library_size = genome.mapped_minus_rDNA)]

# get total length of all exons for each gene
exons_table <-
  width(exons_gr) %>%
  sum(.) %>%
  tibble(gene_id = names(.), width = .) %>% 
  as.data.table(.)

# FPKM
fpkm_long_tb <-
  assay(se) %>%
  as.data.table(., keep.rownames = "gene_id") %>% 
  melt(., 
       id.vars = c("gene_id"),
       variable.name = "sample_id", 
       value.name = "counts") %>% 
  .[, sample_id := str_remove(sample_id, ".genome.Aligned.sortedByCoord.out.bam|.total.bam|.bam")] %>% 
  .[sample_table[, -"bam_path"], on = "sample_id", `:=`(fpm = (counts / round(library_size / 1E6, 6)), 
                                                        stage = i.stage)] %>% 
  .[exons_table, on = "gene_id", fpkm := (fpm / round(width / 1E3, 3))] %>% 
  .[]

# wide and save
fpkm_wide_tb <- 
  fpkm_long_tb %>% 
  dcast(., gene_id ~ sample_id, value.var = "fpkm") %T>% 
  write_csv(., file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".FPKM.csv")))

# average FPKM across stages
fpkm_avg_tb <- 
  fpkm_long_tb %>% 
  .[, list(avg_fpkm = round(mean(fpkm), 3)), by = .(gene_id, stage)] %>%
  dcast(., gene_id ~ stage, value.var = "avg_fpkm") %T>% 
  write_csv(., file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".avgFPKM.csv")))

# long FPKM with mean, SD and standard error for each gene
fpkm_stats <- 
  fpkm_long_tb %>% 
  .[, list(avg_fpkm = mean(fpkm), 
           SD = sd(fpkm), 
           SE = standard.error(fpkm)), 
    by = .(gene_id, stage)] %T>% 
  write_csv(., file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".FPKM_statistics.csv")))

