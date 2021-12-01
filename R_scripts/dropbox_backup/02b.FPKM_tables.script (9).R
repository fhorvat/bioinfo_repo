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
single_end <- as.logical(args$single_end)
threads <- as.numeric(args$threads)
mapped_path <- args$mapped_path
documentation_path <- args$documentation_path
feature_coordinates <- args$feature_coordinates
feature_names <- args$feature_names
grouping_variables <- args$grouping_variables


# counts path
counts_path <- list.files(inpath, ".*\\.multi\\.counts\\.txt$", full.names = T)

# sample table path
sample_table_path <- list.files(documentation_path, ".*sampleTable\\.csv", full.names = T)

# reads stats path
read_stats_path <- file.path(mapped_path, "4_library_size", "library_sizes.txt")

######################################################## READ DATA
# read counts from featureCounts
counts_tb <-
  readr::read_delim(counts_path, delim = "\t", comment = "#") %>%
  set_colnames(., basename(colnames(.)))

# read sample table
sample_table <- data.table::fread(sample_table_path)

# read and clean stats
reads_stats <-
  read_delim(read_stats_path, delim = "\t", col_names = c("sample_id", "library_size")) %>%
  dplyr::filter(!str_detect(sample_id, "21to23nt|24to31nt")) %>%
  as.data.table(.)

######################################################## MAIN CODE
### prepare tables
# get feature coordinates
features_tb <- 
  counts_tb %>% 
  dplyr::select(gene_id = Geneid, seqnames = Chr, start = Start, end = End, width = Length) %>% 
  as.data.table(.)

# counts table
counts_tb %<>%
  dplyr::select(-c(Chr:Length)) %>% 
  dplyr::rename(gene_id = Geneid) %>% 
  as.data.table(.)

# join sample table with stats and tracks
sample_table[reads_stats, on = "sample_id", `:=`(library_size = library_size)]


### calculate FPKMs
# long FPKM/FPM/count table
fpkm_long_tb <-
  counts_tb %>% 
  melt(.,
       id.vars = c("gene_id"),
       variable.name = "sample_id",
       value.name = "counts") %>%
  .[, sample_id := str_remove_all(sample_id, "\\.21to23nt|\\.24to31nt|\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam|\\.total\\.bam|\\.bam")] %>%
  .[sample_table, on = "sample_id"] %>%
  .[features_tb, on = "gene_id", `:=`(width = width)] %>%
  .[, fpm := (counts / round(library_size / 1E6, 6))] %>%
  .[, fpkm := (fpm / round(width / 1E3, 3))] %>%
  .[, c("gene_id", "sample_id", "counts", "fpm", "fpkm", grouping_variables), with = F]

# save long FPKM table
fpkm_long_tb %>%
  .[, -grouping_variables, with = F] %T>%
  readr::write_csv(., file.path(outpath, str_c(feature_names, ".FPKM_long.csv")))

# save wide FPKM table
fpkm_wide_tb <-
  fpkm_long_tb %>%
  dcast(., gene_id ~ sample_id, value.var = "fpkm") %T>%
  write_csv(., file.path(outpath, str_c(feature_names, ".FPKM.csv")))

# average FPKM across stages and/or genotype
fpkm_mean_tb <-
  fpkm_long_tb %>%
  .[, .(fpkm_mean = round(mean(fpkm), 3)), by = c("gene_id", grouping_variables)] %>%
  .[, grouped_variables := do.call(str_c, c(.SD, sep = "_")), .SDcols = grouping_variables] %>%
  dcast(., gene_id ~ grouped_variables, value.var = "fpkm_mean") %T>%
  write_csv(., file.path(outpath, str_c(feature_names, ".FPKM_mean.csv")))

# long FPKM with mean, SD and standard error for each gene
fpkm_stats <-
  fpkm_long_tb %>%
  .[, list(avg_fpkm = mean(fpkm),
           SD = sd(fpkm),
           SE = standard.error(fpkm)),
    by = c("gene_id", grouping_variables)] %T>%
  write_csv(., file.path(outpath, str_c(feature_names, ".FPKM_stats.csv")))
