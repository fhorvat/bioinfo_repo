### INFO: find expression of mature miRNAs in oocyte
### DATE: Fri Sep 13 12:33:04 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/miRNA_expression.miRBase")

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
library(Rsamtools)
library(BiocParallel)
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# miRBase path
mirbase_path <- file.path(genome_dir, "miRBase.22.mm10.20181605.gff3")

# Yang, Garcia-Lopez and Tam small RNA-seq paths
experiment_paths <- c("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/Yang_2016_SciAdv_GSE83581", 
                      "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/GarciaLopez_2015_RNA_GSE59254", 
                      "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/Tam_2008_Nature_GSE10364")

# get list of 21-23nt perfect reads bam files
bam_list <- 
  list.files(str_c(experiment_paths, "/Data/Mapped/STAR_mm10.21_to_23nt.perfect"), pattern = ".*\\.21to23nt\\.bam$", full.names = T) %>% 
  .[str_detect(., "oocyte|MII")]

# get list of library sizes for bam files
library_size_list <- list.files(str_c(experiment_paths, "/Data/Mapped/STAR_mm10.21_to_23nt.perfect/4_library_size"), pattern = "library_sizes\\.txt", full.names = T)

######################################################## READ DATA
# read miRBase
mir_gr <- rtracklayer::import.gff(con = mirbase_path) 

# read and tidy library sizes
library_size_tb <- purrr::map(library_size_list, function(path){
  
  # read, set column names, add experiment
  library_size <- 
    readr::read_delim(path, delim = "\t", col_names = c("bam_name", "library_size")) %>% 
    mutate(experiment = str_extract(path, "Yang|Tam|GarciaLopez"), 
           bam_name = str_c(bam_name, ".bam"))
  
}) %>% 
  bind_rows(.)

######################################################## MAIN CODE
### prepare data
# filter miRBase to include only mature miRNA coordinates
mir_gr <- mir_gr[mcols(mir_gr)$type == "miRNA"]
mcols(mir_gr)$mir_uniq <- str_c("mir", 1:length(mir_gr))
names(mir_gr) <- mcols(mir_gr)$mir_uniq

# get miRBase as table
mir_tb <-
  mir_gr %>% 
  as_tibble() %>% 
  tidyr::unite("coordinates", c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::select(mir_uniq, mir_name = Name, coordinates, strand, width, mir_ID = ID, derives_from = Derives_from)
  
# create sample table
sample_tb <- 
  tibble(bam_path = bam_list) %>% 
  mutate(experiment = str_extract(bam_path, "Yang|Tam|GarciaLopez"), 
         bam_name = basename(bam_path), 
         sample_name = bam_name %>% str_remove(., "\\.SE\\.21to23nt\\.bam") %>% str_c(., experiment, sep = ".")) %>% 
  left_join(., library_size_tb, by = c("bam_name", "experiment"))


### get FPM values of reads over mature miRNAs
# summarizeOverlaps
bamfiles <- Rsamtools::BamFileList(bam_list, yieldSize = 2000000)
se <- GenomicAlignments::summarizeOverlaps(features = mir_gr,
                                           reads = bamfiles,
                                           mode = "IntersectionStrict",
                                           singleEnd = TRUE,
                                           ignore.strand = FALSE)

# get raw counts
mirna_counts <- 
  se %>% 
  assay(.) %>% 
  as_tibble(., rownames = "mir_uniq")

# get FPM normalized values 
mirna_fpm <- 
  mirna_counts %>% 
  tidyr::gather(key = "bam_name", value = "count", -mir_uniq) %>% 
  dplyr::left_join(., sample_tb %>% dplyr::select(bam_name, sample_name, library_size), by = "bam_name") %>% 
  dplyr::mutate(FPM = (count / round(library_size / 1E6, 6))) %>% 
  dplyr::select(mir_uniq, sample_name, FPM) %>% 
  tidyr::spread(key = sample_name, value = FPM) %>% 
  dplyr::left_join(., mir_tb, by = "mir_uniq")

### save as xlsx
# open workbook
wb <- openxlsx::createWorkbook()

# add sheet and write data to the sheet
openxlsx::addWorksheet(wb = wb, sheetName = "oocyte_FPM")
openxlsx::writeData(wb = wb, sheet = "oocyte_FPM", x = mirna_fpm)

# write table to disk
openxlsx::saveWorkbook(wb = wb, 
                       file = file.path(outpath, "mature_miRNA.21to23nt_perfect.FPM_oocyte.miRBase.22.20190913.xlsx"), 
                       overwrite = TRUE)




