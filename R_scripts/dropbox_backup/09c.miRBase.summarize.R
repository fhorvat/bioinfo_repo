### INFO: counts reads over mature miRNAs (from miRbase .gff)
### DATE: Thu Aug 16 17:26:28 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/miRNA_expression")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(xlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# miRbase path
mirbase_path <- list.files(genome_path, pattern = "miRBase.*gff3", full.names = T)

######################################################## READ DATA
# read miRbase gtf
mirbase_gr <- rtracklayer::import.gff(con = mirbase_path) 

######################################################## MAIN CODE
# get ranges of mature miRNA
mirna_gr <- mirbase_gr[mcols(mirbase_gr)$type == "miRNA"]
mcols(mirna_gr) <- mcols(mirna_gr)[, c("Name", "Derives_from")]
mcols(mirna_gr)$unique_name <- str_c(mcols(mirna_gr)$Name, ".", mcols(mirna_gr)$Derives_from)

# create data.frame with miRNA coordinates
mirna_df <- 
  mirna_gr %>% 
  as.data.frame(.) %>% 
  as.tibble(.) %>% 
  tidyr::unite(coordinates, seqnames, start, end, sep = " ") %>% 
  dplyr::select(unique_name, coordinates, strand) %>% 
  dplyr::filter(!duplicated(unique_name))

### read summarizedOverlaps .RDS for all experiments and writes results table
# experiment paths
experiment_paths <- 
  "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/NIH3T3_transfected.2018/Data/Mapped/STAR_mm10.all_plasmid" %>% 
  magrittr::set_names(., str_remove(., "Data.*$") %>% basename(.))

# list of experiments
experiment_list <- names(experiment_paths)

# loop through experiments
experiment <- experiment_list[1]

# summarizedOverlaps .RDS path
se_path <- file.path(outpath, str_c("miRbase.", experiment, ".se.RDS"))

# library size path
library_size_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/library_size", str_c("library_size.counts.", experiment, ".csv"))

# read summarizeOverlaps .RDS
se <- readRDS(file = se_path)

# read library size table
library_size <- 
  readr::read_csv(library_size_path) %>% 
  tidyr::gather(mismatch, library_size, -sample_id) %>% 
  tidyr::separate(mismatch, c("mismatch", "read_group"), sep = "\\.") %>% 
  tidyr::unite(sample_id, sample_id, mismatch, sep = ".") %>% 
  tidyr::spread(read_group, library_size) %>% 
  dplyr::select(sample_id, library_size = `21to23nt_reads`)

# clean miRNA summarizeOverlaps
mirna_clean <- 
  assay(se) %>% 
  as.tibble(.) %>% 
  magrittr::set_colnames(str_remove(colnames(.), ".bam")) %>% 
  dplyr::mutate(mirna_id = mcols(mirna_gr)$unique_name) %>% 
  dplyr::select(mirna_id, everything()) %>% 
  dplyr::filter(!duplicated(mirna_id)) %>% 
  tidyr::gather(sample_id, count, -mirna_id) %>% 
  tidyr::separate(sample_id, c("sample_id", "mismatch"), sep = ".mis_") %>% 
  dplyr::mutate(mismatch = str_c("mis_", mismatch)) %>% 
  tidyr::spread(mismatch, count) %>% 
  dplyr::mutate_at(vars(starts_with("mis")), funs(replace(., is.na(.), 0))) %>% 
  dplyr::mutate(mis_1 = mis_1 + mis_0, 
                mis_2 = mis_2 + mis_1) %>% 
  dplyr::arrange(sample_id) %>% 
  tidyr::gather(mismatch_allowed, count, mis_0:mis_2) 

# counts
mirna_counts <- 
  mirna_clean %>%
  tidyr::unite(tmp, mismatch_allowed, sample_id, sep = ".") %>% 
  tidyr::spread(tmp, count) %>%
  dplyr::mutate_at(vars(starts_with("mis_")), funs(replace(., is.na(.), 0))) %>% 
  dplyr::left_join(., mirna_df %>% dplyr::select(-strand), by = c("mirna_id" = "unique_name")) %>% 
  dplyr::select(mirna_id, coordinates, everything())

# RPM
mirna_rpm <- 
  mirna_clean %>%
  tidyr::unite(sample_id, sample_id, mismatch_allowed, sep = ".") %>% 
  dplyr::left_join(library_size, by = "sample_id") %>% 
  dplyr::mutate(library_size = library_size / 1e6, 
                rpm = round((count / library_size), 3), 
                tmp = str_c(str_extract(sample_id, "mis_[0-2]{1}"), ".",  
                            str_remove(sample_id, "\\.mis_[0-2]{1}"))) %>% 
  dplyr::select(mirna_id, tmp, rpm) %>%
  tidyr::spread(tmp, rpm) %>%
  dplyr::mutate_at(vars(starts_with("mis_")), funs(replace(., is.na(.), 0))) %>% 
  dplyr::left_join(., mirna_df %>% dplyr::select(-strand), by = c("mirna_id" = "unique_name")) %>% 
  dplyr::select(mirna_id, coordinates, everything())

# write RPMs to separate sheets in xlsx 
write.xlsx(x = mirna_rpm %>% 
             dplyr::select_at(vars("mirna_id", "coordinates", starts_with("mis_0"))) %>% 
             magrittr::set_colnames(., str_remove(colnames(.), "mis_[0-2]{1}_")) %>% 
             dplyr::arrange(desc(rowSums(.[, -c(1, 2)]))) %>% 
             as.data.frame(.), 
           file = file.path(outpath, str_c("mature_miRNA.", experiment, ".21to23nt.RPM.xlsx")), 
           sheetName = "allowed_0_mismatches", 
           row.names = FALSE)

write.xlsx(x = mirna_rpm %>% 
             dplyr::select_at(vars("mirna_id", "coordinates", starts_with("mis_1"))) %>% 
             magrittr::set_colnames(., str_remove(colnames(.), "mis_[0-2]{1}_")) %>% 
             dplyr::arrange(desc(rowSums(.[, -c(1, 2)]))) %>%
             as.data.frame(.), 
           file = file.path(outpath, str_c("mature_miRNA.", experiment, ".21to23nt.RPM.xlsx")), 
           sheetName = "allowed_1_mismatches", 
           append = TRUE, 
           row.names = FALSE)

write.xlsx(x = mirna_rpm %>% 
             dplyr::select_at(vars("mirna_id", "coordinates", starts_with("mis_2"))) %>% 
             magrittr::set_colnames(., str_remove(colnames(.), "mis_[0-2]{1}_")) %>% 
             dplyr::arrange(desc(rowSums(.[, -c(1, 2)]))) %>% 
             as.data.frame(.), 
           file = file.path(outpath, str_c("mature_miRNA.", experiment, ".21to23nt.RPM.xlsx")), 
           sheetName = "allowed_2_mismatches", 
           append = TRUE, 
           row.names = FALSE)

# write counts to separate sheets in xlsx 
write.xlsx(x = mirna_counts %>% 
             dplyr::select_at(vars("mirna_id", "coordinates", starts_with("mis_0"))) %>% 
             magrittr::set_colnames(., str_remove(colnames(.), "mis_[0-2]{1}_")) %>% 
             dplyr::arrange(desc(rowSums(.[, -c(1, 2)]))) %>% 
             as.data.frame(.), 
           file = file.path(outpath, str_c("mature_miRNA.", experiment, ".counts.xlsx")), 
           sheetName = "allowed_0_mismatches", 
           row.names = FALSE)

write.xlsx(x = mirna_counts %>% 
             dplyr::select_at(vars("mirna_id", "coordinates", starts_with("mis_1"))) %>% 
             magrittr::set_colnames(., str_remove(colnames(.), "mis_[0-2]{1}_")) %>% 
             dplyr::arrange(desc(rowSums(.[, -c(1, 2)]))) %>%
             as.data.frame(.), 
           file = file.path(outpath, str_c("mature_miRNA.", experiment, ".counts.xlsx")), 
           sheetName = "allowed_1_mismatches", 
           append = TRUE, 
           row.names = FALSE)

write.xlsx(x = mirna_counts %>% 
             dplyr::select_at(vars("mirna_id", "coordinates", starts_with("mis_2"))) %>% 
             magrittr::set_colnames(., str_remove(colnames(.), "mis_[0-2]{1}_")) %>% 
             dplyr::arrange(desc(rowSums(.[, -c(1, 2)]))) %>% 
             as.data.frame(.), 
           file = file.path(outpath, str_c("mature_miRNA.", experiment, ".counts.xlsx")), 
           sheetName = "allowed_2_mismatches", 
           append = TRUE, 
           row.names = FALSE)
