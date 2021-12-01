#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO:
### DATE: Thu Apr 25 16:16:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd(".")

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
# calculate standard error of the mean FPKM value
standard.error <- function(x) {
  sqrt(var(x) / length(x))
}

######################################################## PATH VARIABLES
### in and out
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# arguments from command line
experiment <- args$experiment
ensembl_version <- args$ensembl_version
genome_path <- args$genome_path
mapped_path <- args$mapped_path
threads <- args$threads
grouping_variables <- args$grouping_variables

ensembl_version="93"
genome_path="/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"
mapped_path="STAR_mesAur1"
threads="1"
grouping_variables="genotype"
experiment="hamster_oocyte_Mov10l"
mapped_path="/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq/Data/Mapped/STAR_mesAur1"

### other experiment paths
# set base experiment path
base_path <- file.path(mapped_path, "../..")

# documentation path
documentation_path <- file.path(base_path, "Documentation")

# analysis path
analysis_path <- inpath


### documentation
# sample table path
sample_table_path <- list.files(documentation_path, ".*sampleTable\\.csv", full.names = T)

# stats and tracks path
stats_and_tracks_path <- list.files(mapped_path, ".*\\.stats_and_tracks\\.csv", full.names = T)

# summarizedExperiment path
se_path <- list.files(path = analysis_path, pattern = ".*\\.se\\.RDS$", full.names = T)


### genome
# reduced exons path
exons_path <- file.path(inpath, "joined_coverage.perfect_multimappers.coordinates.csv")

######################################################## READ DATA
# read sample table
sample_table <- data.table::fread(sample_table_path)

# read stats and tracks table
stats_and_tracks <- data.table::fread(stats_and_tracks_path)

# summarizedExperiment
se <- readRDS(se_path)

# read ENSEMBL reduced exons
exons_gr <- 
  readr::read_csv(file = exons_path) %>% 
  GRanges(.)
names(exons_gr) <- mcols(exons_gr)$gene_id
exons_gr <- split(exons_gr, names(exons_gr))

######################################################## MAIN CODE
# join sample table with stats and tracks
sample_table[stats_and_tracks, on = "sample_id", `:=`(library_size = genome.mapped_minus_rDNA)]

# get total length of all exons for each gene
exons_table <-
  width(exons_gr) %>%
  sum(.) %>%
  tibble(gene_id = names(.), width = .) %>%
  as.data.table(.)


### calculate FPKMs
# long FPKM/FPM/count table
fpkm_long_tb <-
  assay(se) %>%
  as.data.table(., keep.rownames = "gene_id") %>%
  melt(.,
       id.vars = c("gene_id"),
       variable.name = "sample_id",
       value.name = "counts") %>%
  .[, sample_id := str_remove(sample_id, "\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam|\\.total\\.bam|\\.bam")] %>%
  .[sample_table, on = "sample_id"] %>%
  .[exons_table, on = "gene_id"] %>%
  .[, fpm := (counts / round(library_size / 1E6, 6))] %>%
  .[, fpkm := (fpm / round(width / 1E3, 3))] %>%
  .[, c("gene_id", "sample_id", "counts", "fpm", "fpkm", grouping_variables), with = F]

# save long FPKM table
fpkm_long_tb %>%
  .[, -grouping_variables, with = F] %>%
  readr::write_csv(., file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".FPKM_long.csv")))

# save wide FPKM table
fpkm_wide_tb <-
  fpkm_long_tb %>%
  dcast(., gene_id ~ sample_id, value.var = "fpkm") %T>%
  write_csv(., file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".FPKM.csv")))

# average FPKM across stages and/or genotype
fpkm_mean_tb <-
  fpkm_long_tb %>%
  .[, .(fpkm_mean = round(mean(fpkm), 3)), by = c("gene_id", grouping_variables)] %>%
  .[, grouped_variables := do.call(str_c, c(.SD, sep = "_")), .SDcols = grouping_variables] %>%
  dcast(., gene_id ~ grouped_variables, value.var = "fpkm_mean") %>%
  write_csv(., file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".FPKM_mean.csv")))

# long FPKM with mean, SD and standard error for each gene
fpkm_stats <-
  fpkm_long_tb %>%
  .[, list(avg_fpkm = mean(fpkm),
           SD = sd(fpkm),
           SE = standard.error(fpkm)),
    by = c("gene_id", grouping_variables)] %T>%
  write_csv(., file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".FPKM_stats.csv")))
