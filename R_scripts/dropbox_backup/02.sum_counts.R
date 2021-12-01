### INFO: read smallRNA seq bam file, get counts over exons
### DATE: Wed May 16 02:54:16 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# mouse
animal <- "mouse"
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/mouse")
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# # cow
# animal <- "cow"
# setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/cow")
# genome_path <- "/common/DB/genome_reference/cow/bosTau8.UMD3.1.GCA_000003055.4"

# # pig
# animal <- "pig"
# setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/pig")
# genome_path <- "/common/DB/genome_reference/pig/susScr11.Sscrofa11.1.GCA_000003025.6"

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(purrr)
library(data.table)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)

######################################################## SOURCE FILES
# lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
# source(file.path(lib_path, "wideScreen.R"))
# source(file.path(lib_path, "headt.R"))
# source(file.path(lib_path, "asdf.R"))
# source(file.path(lib_path, "bamToGRangesList.R"))
# wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main outpath
outpath <- getwd()

# set main inpath
inpath <- getwd()

# tables list
table_list <- list.files(inpath, pattern = "smallRNA.perfect.21to23.summary.counts.*csv", full.names = T)

# get paths of reduced exons
exons_path <- list.files(inpath, "ensembl.91.*filt.rmsk_miRNA.RDS", full.names = T)

# gene info path
info_path <- list.files(genome_path, "ensembl.91.*[0-9]{6}.UCSCseqnames.geneInfo.csv", full.names = T)

# library size path
libsize_list <-
  list.files(file.path(inpath, "../data_sets", animal), pattern = ".*lib_size.hist.txt",
             recursive = T, full.names = T) %>%
  tibble(lib_path = .) %>%
  dplyr::mutate(sample = basename(lib_path) %>% str_remove_all(., ".SE.lib_size.hist.txt|.SE_cl.lib_size.hist.txt|.SE.*"),
                experiment = str_remove(lib_path, "\\/Links.*") %>% basename(.) %>% str_remove(., "(?<=[0-9])_.*"))

######################################################## READ DATA
# read info about genes
genes_info <- readr::read_csv(info_path)

# read library sizes to one data.frame
libsize_df <-
  purrr::pmap(libsize_list, function(lib_path, sample, experiment) {
    
    # read library data.frame from files defined in libsize_list data.frame
    lib_df <-
      data.table::fread(lib_path) %>%
      as.tibble(.) %>%
      magrittr::set_colnames(c("lib_size", "seq_length")) %>%
      dplyr::mutate(ID = str_c(experiment, sample, sep = "."))
    
  }) %>%
  dplyr::bind_rows(.) %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(lib_size_all = sum(lib_size),
                   lib_size_19_24 = sum(lib_size[seq_length >= 19 & seq_length <= 24]),
                   lib_size_21_23 = sum(lib_size[seq_length >= 21 & seq_length <= 23])) %>%
  dplyr::mutate_at(vars(starts_with("lib_size")), funs(round(. / 1E6, 5)))

######################################################## MAIN CODE
# create empty table with all categories in order
category_df <-
  tibble(read_group = rep(c("rRNA", "SINE", "LINE", "LTR", "other_repeat", "miRNA", "lincRNA", "protein_coding", "annotated_pseudogene", "other"), each = 2),
         strand = rep(c("sense", "antisense"), length(read_group) / 2)) %>%
  dplyr::bind_rows(., tibble(read_group = "not_annotated", strand = NA))

# read all tables, join with all categories
count_df <-
  purrr::map(table_list, function(table_in){
    
    readr::read_csv(table_in) %>%
      dplyr::full_join(., category_df, by = c("read_group", "strand")) %>%
      dplyr::mutate(sample = replace(sample, is.na(sample), sample[1]),
                    experiment = replace(experiment, is.na(experiment), experiment[1]),
                    count = replace(count, is.na(count), 0))
    
  }) %>%
  dplyr::bind_rows(.) %>%
  dplyr::filter((experiment != "Tam_2008" | sample == "s_oocyte_19to24")) %>%
  tidyr::unite(col = "ID", experiment, sample, sep = ".") %>%
  dplyr::left_join(., libsize_df, by = "ID") %>%
  dplyr::mutate(fpm_all = count / lib_size_all,
                fpm_19_24 = count / lib_size_19_24,
                fpm_21_23 = count / lib_size_21_23) %>%
  dplyr::mutate_at(vars(starts_with("fpm")), funs(round(., 3)))

# calculate average count/FPMs
count_avg <-
  count_df %>%
  dplyr::group_by(read_group, strand) %>%
  dplyr::summarise(avg_count = mean(count) %>% round(., 2),
                   avg_fpm_all = mean(fpm_all) %>% round(., 2),
                   avg_fpm_19_24 = mean(fpm_19_24) %>% round(., 2),
                   avg_fpm_21_23 = mean(fpm_21_23) %>% round(., 2))

# gather, unite, spread
count_df %<>%
  dplyr::select_at(vars(-starts_with("lib_size"))) %>%
  tidyr::gather(variable, value, -(c(read_group, strand, ID))) %>%
  tidyr::unite(temp, ID, variable) %>%
  tidyr::spread(temp, value) %>%
  dplyr::left_join(., count_avg, by = c("read_group", "strand")) %>%
  dplyr::right_join(., category_df, by = c("read_group", "strand")) %T>%
  readr::write_csv(x = ., path = file.path(outpath, str_c("smallRNA.perfect.21to23.summary.counts.all.", animal, ".csv")))
