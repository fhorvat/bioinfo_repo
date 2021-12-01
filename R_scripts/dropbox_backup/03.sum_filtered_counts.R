### INFO: read smallRNA seq bam file, get counts over exons
### DATE: Wed May 16 02:54:16 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# # mouse
# animal <- "mouse"
# setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/mouse")
# genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# cow
animal <- "cow"
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/cow")
genome_path <- "/common/DB/genome_reference/cow/bosTau8.UMD3.1.GCA_000003055.4"

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

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "bamToGRangesList.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main outpath
outpath <- getwd()

# set main inpath
inpath <- getwd()

# tables list 
table_list <- list.files(inpath, pattern = "smallRNA.perfect.21to23.counts.*csv", full.names = T)

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
# read all tables
count_df <- 
  purrr::map(table_list, readr::read_csv) %>% 
  dplyr::bind_rows(.)

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
  dplyr::bind_rows(.)

######################################################## MAIN CODE
# filter library sizes
libsize_df %<>% 
  dplyr::filter(seq_length >= 19, seq_length <= 24) %>% 
  dplyr::group_by(ID) %>% 
  dplyr::summarise(lib_size = sum(lib_size)) %>% 
  dplyr::mutate(lib_size_million = round(lib_size / 1E6, 6))

# get coordinates and info about features
coord_df <- 
  genes_info %>% 
  dplyr::mutate(coordinates = str_c(seqnames, start, end, sep = " ")) %>% 
  dplyr::select(gene_id, gene_name, coordinates, strand, gene_biotype, gene_description)

### get results
# filter count tables
count_df %<>% 
  dplyr::filter((experiment != "Tam_2008" | sample == "s_oocyte_19to24")) %>% 
  tidyr::unite(col = "ID", experiment, sample, sep = ".") %>% 
  dplyr::left_join(., libsize_df %>% dplyr::select(ID, lib_size_million), by = "ID") %>% 
  dplyr::mutate(fpm_sense = count_sense / lib_size_million, 
                fpm_antisense = count_antisense / lib_size_million)

# calculate average count per exon per aligned strand
count_avg <- 
  count_df %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::mutate(count_all = count_sense + count_antisense, 
                fpm_all = fpm_sense + fpm_antisense) %>% 
  dplyr::summarise(avg_count.all = mean(count_all) %>% round(., 2),
                   avg_count.sense = mean(count_sense) %>% round(., 2), 
                   avg_count.antisense = mean(count_antisense) %>% round(., 2),
                   avg_fpm.all = mean(fpm_all) %>% round(., 2),
                   avg_fpm.sense = mean(fpm_sense) %>% round(., 2), 
                   avg_fpm.antisense = mean(fpm_antisense) %>% round(., 2))

# gather, unite, spread
count_df %<>% 
  dplyr::select(-lib_size_million) %>% 
  tidyr::gather(key = ID:gene_id, value = count, count_sense, count_antisense, fpm_sense, fpm_antisense) %>% 
  dplyr::rename(alignment_sense = `ID:gene_id`) %>% 
  tidyr::unite(col = ID, ID, alignment_sense, sep = ".") %>% 
  tidyr::spread(key = ID, value = count) %>%
  dplyr::left_join(., count_avg, by = "gene_id")

# split to raw and normalized tables
count_df_out <-
  count_df %>% 
  dplyr::select_at(.vars = vars(gene_id, contains("count"))) %>% 
  dplyr::left_join(., coord_df, by = "gene_id") %>% 
  dplyr::select(gene_id, gene_name, coordinates, strand, 
                avg_count.all, avg_count.sense, avg_count.antisense, 
                everything()) %>% 
  dplyr::arrange(desc(avg_count.all)) %>% 
  dplyr::filter(avg_count.all > 0) %T>% 
  readr::write_csv(x = ., path = file.path(outpath, str_c("smallRNA.perfect.21to23.all.raw.", animal, ".csv")))

# split to raw and normalized tables
fpm_df_out <-
  count_df %>% 
  dplyr::select_at(.vars = vars(gene_id, contains("fpm"))) %>% 
  dplyr::left_join(., coord_df, by = "gene_id") %>% 
  dplyr::select(gene_id, gene_name, coordinates, strand, 
                avg_fpm.all, avg_fpm.sense, avg_fpm.antisense, 
                everything()) %>% 
  dplyr::arrange(desc(avg_fpm.all)) %>% 
  dplyr::filter(avg_fpm.all > 0) %T>% 
  readr::write_csv(x = ., path = file.path(outpath, str_c("smallRNA.perfect.21to23.all.fpm.", animal, ".csv")))



