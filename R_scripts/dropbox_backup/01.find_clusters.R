### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters")

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

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"

# documentation path
documentation_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq/Documentation")

# sample table path
sample_tb_path <- list.files(documentation_path, ".*sampleTable.csv$")

# mapped path
mapped_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq/Data/Mapped/STAR_mesAur1.new/8_merged_replicates")

# read counts path
reads_tb_path <- list.files(file.path(mapped_path, "../4_library_size"), "library_sizes.txt", full.names = T, recursive = T)

# list piRNA .bams
bams_path <- 
  list.files(mapped_path, pattern = "s_testis_.*13dpp.*24to31nt\\.bam$", full.names = T) %>% 
  .[!str_detect(., "HET")]

######################################################## READ DATA
# read sample table
sample_tb <- readr::read_csv(sample_tb_path)

# read read counts
reads_tb <- readr::read_delim(reads_tb_path, delim = "\t", col_names = c("sample_id", "library_size"))

# read bams
bams_list <- purrr::map(bams_path, function(path){
  
  # read bam as GenomicRanges
  bam_gr <- GenomicAlignments::readGAlignments(path)
  
}) %>% 
  set_names(., basename(bams_path) %>% str_remove(., "\\.bam$"))

######################################################## MAIN CODE
# clean and filter read stats
reads_tb_tidy <- 
  reads_tb %>% 
  dplyr::filter(str_detect(sample_id, "19to32nt"), 
                str_detect(sample_id, "13dpp"), 
                !str_detect(sample_id, "HET")) %>% 
  dplyr::mutate(sample_id = str_extract(sample_id, "KO|WT")) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::summarise(library_size = sum(library_size)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(library_size = (library_size / 10e5))

# read coordinates in a table
reads_tb_list <- purrr::map(names(bams_list), function(name){
  
  # get table, find starts of piRNAs, count
  reads_tb <- 
    bams_list[[name]] %>% 
    as_tibble(.) %>% 
    dplyr::mutate(coordinates = ifelse(strand == "+", 
                                       str_c(seqnames, start, sep = ":"), 
                                       str_c(seqnames, end, sep = ":"))) %>% 
    dplyr::group_by(coordinates) %>% 
    dplyr::summarise(count = n()) %>% 
    dplyr::arrange(desc(count))
  
  # return
  return(reads_tb)
  
}) %>% 
  set_names(., names(bams_list))

# get RPK values of piRNAs
reads_tb_rpm <- 
  reads_tb_list[["s_testis_Mov10l_WT_13dpp.24to31nt"]] %>% 
  dplyr::filter(count >= 2) %>% 
  dplyr::left_join(., reads_tb_list[["s_testis_Mov10l_KO_13dpp.24to31nt"]], by = "coordinates") %>% 
  dplyr::rename(count_WT = count.x, count_KO = count.y) %>% 
  dplyr::mutate(count_KO = replace(count_KO, is.na(count_KO), 0)) %>% 
  tidyr::pivot_longer(cols = -coordinates, names_prefix = "count_", names_to = "sample_id", values_to = "count") %>% 
  dplyr::left_join(., reads_tb_tidy, by = "sample_id") %>% 
  dplyr::mutate(rpm = round((count / library_size), 3)) %>% 
  tidyr::pivot_wider(id_cols = coordinates, names_from = "sample_id", values_from = "rpm", names_prefix = "rpm_")

# get piRNAs which are downregulated in KO
reads_tb_filt <- 
  reads_tb_rpm %>% 
  dplyr::filter(rpm_KO <= (rpm_WT / 4)) %>% 
  dplyr::filter(!str_detect(coordinates, "^MT:"))

# save
readr::write_csv(reads_tb_filt, file.path(outpath, "piRNA_reads.individual_reads.13dpp_Mov10l1_KO.downregulated.csv"))


### get full coordinates
# get coordinates, overlap with actual piRNAs, get full coordinates
reads_gr <- 
  reads_tb_filt %>% 
  tidyr::separate(coordinates, into = c("seqnames", "start"), sep = ":") %>% 
  dplyr::mutate(end = start) %>% 
  GRanges(.) 

# find overlaps
overlaps <- findOverlaps(reads_gr, bams_list[["s_testis_Mov10l_WT_13dpp.24to31nt"]], ignore.strand = T)

# get hits
pirna_gr <- 
  bams_list[["s_testis_Mov10l_WT_13dpp.24to31nt"]] %>% 
  GRanges(.) %>% 
  .[subjectHits(overlaps)] %>% 
  unique(.) %>% 
  GenomicRanges::reduce(., ignore.strand = T)

# count overlaps between hits and reads
wt_counts <- GenomicRanges::countOverlaps(pirna_gr, bams_list[["s_testis_Mov10l_WT_13dpp.24to31nt"]], ignore.strand = T, minoverlap = 20)
ko_counts <- GenomicRanges::countOverlaps(pirna_gr, bams_list[["s_testis_Mov10l_KO_13dpp.24to31nt"]], ignore.strand = T, minoverlap = 20)

# create table, join with counts, calculate RPM 
pirna_tb <- 
  pirna_gr %>% 
  as_tibble(.) %>% 
  dplyr::select(-c("width", "strand")) %>% 
  tidyr::unite(col = "coordinates", seqnames:end, sep = " ") %>% 
  dplyr::mutate(count_WT = wt_counts, 
                count_KO = ko_counts) %>% 
  tidyr::pivot_longer(cols = -coordinates, names_prefix = "count_", names_to = "sample_id", values_to = "count") %>% 
  dplyr::left_join(., reads_tb_tidy, by = "sample_id") %>% 
  dplyr::mutate(rpm = round((count / library_size), 3)) %>% 
  tidyr::pivot_wider(id_cols = coordinates, names_from = "sample_id", values_from = "rpm", names_prefix = "rpm_") %>% 
  dplyr::filter(rpm_KO <= (rpm_WT / 4)) %>% 
  arrange(-rpm_WT)

# save 
readr::write_csv(pirna_tb, file.path(outpath, "piRNA_reads.joined_reads.13dpp_Mov10l1_KO.downregulated.csv"))
  
  
  