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
features_coordinates <- args$features_coordinates
features_name <- args$features_name
genes_info_path <- args$genes_info_path
grouping_variables <- args$grouping_variables
counts_path <- args$counts_path

# sample table path
sample_table_path <- list.files(documentation_path, ".*sampleTable\\.csv", full.names = T)

# reads stats path
reads_stats_path <- file.path(mapped_path, "../4_library_size/library_sizes.txt")

######################################################## READ DATA
# read counts from featureCounts
counts_tb <-
  readr::read_delim(counts_path, delim = "\t", comment = "#") %>%
  set_colnames(., basename(colnames(.))) %>%
  dplyr::mutate(Geneid = make.unique(Geneid))

# read sample table
sample_table <- data.table::fread(sample_table_path)

# read stats and tracks table
reads_stats <-
  readr::read_delim(reads_stats_path, delim = "\t", col_names = c("sample_id", "library_size")) %>%
  dplyr::filter(!is.na(sample_id),
                str_detect(sample_id, "19to32nt$")) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt$")) %>%
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
fpm_long_tb <-
  counts_tb %>%
  melt(.,
       id.vars = c("gene_id"),
       variable.name = "sample_id",
       value.name = "counts") %>%
  .[, sample_id := str_remove_all(sample_id, "\\.21to23nt|\\.24to31nt|\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam|\\.total\\.bam|\\.bam|\\.mis_0")] %>%
  .[sample_table, on = "sample_id"] %>%
  .[features_tb, on = "gene_id", `:=`(width = width,
                                      coordinates = str_c(seqnames, start, end, sep = " "))] %>%
  .[, fpm := (counts / round(library_size / 1E6, 6))] %>%
  .[, fpkm := (fpm / round(width / 1E3, 3))] %>%
  .[, c("gene_id", "sample_id", "counts", "fpm", "fpkm", grouping_variables, "coordinates"), with = F] %>%
  .[!is.na(gene_id)]

# save long FPM table
fpm_long_tb %>%
  .[, -grouping_variables, with = F] %T>%
  readr::write_csv(., file.path(outpath, str_c(features_name, ".FPM_long.csv")))

# save wide FPM table
fpm_wide_tb <-
  fpm_long_tb %>%
  dcast(., gene_id + coordinates ~ sample_id, value.var = "fpm") %T>%
  write_csv(., file.path(outpath, str_c(features_name, ".FPM.csv")))

# average FPM across stages and/or genotype
fpm_mean_tb <-
  fpm_long_tb %>%
  .[, .(fpm_mean = round(mean(fpm), 3),
        coordinates = unique(coordinates)), by = c("gene_id", grouping_variables)] %>%
  .[, grouped_variables := do.call(str_c, c(.SD, sep = "_")), .SDcols = grouping_variables] %>%
  dcast(., gene_id + coordinates ~ grouped_variables, value.var = "fpm_mean") %T>%
  write_csv(., file.path(outpath, str_c(features_name, ".FPM_mean.csv")))

# long FPM with mean, SD and standard error for each gene
fpm_stats <-
  fpm_long_tb %>%
  .[, list(avg_fpm = mean(fpm),
           SD = sd(fpm),
           SE = standard.error(fpm),
           coordinates = unique(coordinates)),
    by = c("gene_id", grouping_variables)] %T>%
  write_csv(., file.path(outpath, str_c(features_name, ".FPM_stats.csv")))
