#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: summarize read stats
### DATE: 06. 10. 2017.
### AUTHOR: Filip Horvat

######################################################## WORKING DIRECTORY

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tibble)

######################################################## SOURCE FILES
# read STAR log parser
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "read_STAR.Log.final.out.R"))

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# outpath
outpath <- getwd()

######################################################## READ DATA
# # get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# list of read stats, genome logs and merged logs
stats_list <- args$read_stats
genome_logs_list <- args$genome_logs
merged_logs_list <- args$merged_logs

######################################################## MAIN CODE
# sample list
sample_list <-
  stats_list %>%
  basename(.) %>%
  stringr::str_remove(., pattern = "\\.read_stats\\.txt$")

# get all posible combinations of alignments (genome/genome.merged)
alignments_all <-
  expand.grid(sample_list, ".", c("genome", "merged")) %>%
  do.call(stringr::str_c, .) %>%
  unique(.) %>%
  tibble(log_id = .)

# reads class
log_class_df <-
  lapply(stats_list, readr::read_delim, delim = "\t") %>%
  dplyr::bind_rows(.)

# log from STAR
log_star <-
  lapply(c(genome_logs_list, merged_logs_list), read_STAR.Log.final.out, reshape = T) %>%
  dplyr::bind_rows(.) %>%
  dplyr::full_join(., alignments_all, by = "log_id") %>%
  dplyr::mutate_all(.funs = funs(replace(., is.na(.), 0))) %>%
  mutate(alignment = str_extract(log_id, pattern = "genome$|merged$"),
         sample_id = str_remove(log_id, pattern = "\\.genome.*|\\.merged.*")) %>%
  dplyr::select(-log_id) %>%
  dplyr::select(sample_id, alignment, everything()) %>%
  magrittr::set_colnames(stringr::str_replace_all(colnames(.), "_reads|number_of_", ""))

# mapped reads
log_mapped <-
  log_star %>%
  dplyr::select(-c(unmapped, input, alignment)) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise_all(.funs = sum) %>%
  dplyr::mutate(genome.mapped_total = uniquely_mapped + mapped_to_multiple_loci, 
                mapped_total = genome.mapped_total) %>%
  dplyr::rename(genome.uniquely_mapped = uniquely_mapped,
                genome.mapped_to_multiple_loci = mapped_to_multiple_loci)

# unmapped reads
log_unmapped <-
  log_star %>%
  dplyr::select(-c(uniquely_mapped, mapped_to_multiple_loci)) %>%
  tidyr::gather(category, number_of_reads, -c(sample_id, alignment)) %>%
  tidyr::unite(temp, alignment, category, sep = ".") %>%
  tidyr::spread(temp, number_of_reads) %>%
  dplyr::mutate(discarded_while_merging = genome.unmapped - merged.input,
                genome.unmapped_total = discarded_while_merging + merged.unmapped)

# combine mapped and unmapped, combine with class of reads mapped to genome
log_final <-
  dplyr::left_join(log_mapped, log_unmapped, by = "sample_id") %>%
  dplyr::left_join(., log_class_df, by = "sample_id") %>%
  dplyr::mutate(unmapped_total = genome.input - mapped_total) %>%
  dplyr::select(sample_id, raw_input = genome.input, mapped_total, unmapped_total,
                genome.mapped_total, genome.uniquely_mapped, genome.mapped_to_multiple_loci, genome.mapped_minus_rDNA, genome.unmapped_total,
                total, rRNA, repeats, exon, other) %>%
  dplyr::mutate_all(.funs = funs(as.character(.))) %>%
  readr::write_delim(x = ., path = "log.read_stats.txt", delim = "\t")
