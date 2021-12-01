### INFO: counts reads over mature miRNAs (from miRbase .gff)
### DATE: Wed Nov 28 15:03:32 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Analysis")

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

# summarizedOverlaps .RDS path
se_path <- file.path(inpath, "miRbase.Eliska_mESC_MosIR.se.RDS")

# library size path
library_size_path <- file.path(inpath, "library_size.Eliska_mESC_MosIR.csv")

######################################################## READ DATA
# read miRbase gtf
mirbase_gr <- rtracklayer::import.gff(con = mirbase_path) 

# read summarizeOverlaps .RDS
se <- readRDS(file = se_path)

# read library size
library_size <- readr::read_csv(library_size_path)

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

# clean miRNA summarizeOverlaps
mirna_clean <- 
  assay(se) %>% 
  as.tibble(.) %>% 
  magrittr::set_colnames(str_remove(colnames(.), "\\.SE\\.mis_0\\.18to30nt\\.bam")) %>% 
  dplyr::mutate(mirna_id = mcols(mirna_gr)$unique_name) %>% 
  dplyr::select(mirna_id, everything()) %>% 
  dplyr::filter(!duplicated(mirna_id)) %>% 
  tidyr::gather(sample_id, alignment_count, -mirna_id) %>% 
  dplyr::left_join(., mirna_df %>% dplyr::select(-strand), by = c("mirna_id" = "unique_name")) 
  
# counts
mirna_counts <- 
  mirna_clean %>% 
  tidyr::spread(sample_id, alignment_count) %>% 
  dplyr::select(mirna_id, coordinates, everything())

# RPM
mirna_rpm <- 
  mirna_clean %>%
  dplyr::left_join(library_size, by = "sample_id") %>% 
  dplyr::mutate(library_size.mis_0.18to30nt = library_size.mis_0.18to30nt / 1e6, 
                rpm = round(alignment_count / library_size.mis_0.18to30nt, 3)) %>% 
  dplyr::select(mirna_id, sample_id, rpm, coordinates) %>%
  tidyr::spread(sample_id, rpm) %>%
  dplyr::select(mirna_id, coordinates, everything())

### write to separate sheets in xlsx 
# counts
write.xlsx(x = as.data.frame(mirna_counts), 
           file = file.path(outpath, "mature_miRNA.Eliska_mESC_MosIR.expression.xlsx"), 
           sheetName = "allowed_0_mismatches.counts", 
           row.names = FALSE)

# RPM
write.xlsx(x = as.data.frame(mirna_rpm), 
           file = file.path(outpath, "mature_miRNA.Eliska_mESC_MosIR.expression.xlsx"), 
           sheetName = "allowed_0_mismatches.RPMs", 
           row.names = FALSE, 
           append = TRUE)

