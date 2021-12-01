### INFO: read smallRNA seq bam file, get counts over exons
### DATE: Mon May 14 12:21:35 2018
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
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# get paths of reduced exons
exons_path <- list.files(inpath, "ensembl.91.*filt.rmsk_miRNA.RDS", full.names = T)

# get bam file path, name, experiment
# bam_path <- "%BAM_PATH"
# bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/data_sets/cow/Roovers_2015_CellRep_GSE64942/Filtered/s_bosTau_GV_r1.SE.perfect.bam"
bam_name <- basename(bam_path) %>% str_remove(., ".SE.perfect.bam")
experiment_name <- bam_path %>% str_remove(., "/Filtered/.*") %>% basename(.) %>% str_remove(., "(?<=[0-9])_.*")

######################################################## READ DATA
# exons by genes
exons_gr <- readRDS(file = exons_path)

# read bam as ranges
bam_gr <- 
  bamToGRangesList(bam_path) %>% 
  unlist(.)

######################################################## MAIN CODE
### prepare bam ranges to count
# get reads with width between 21 and 23 nt
bam_gr_filt <- bam_gr[(width(bam_gr) >= 21) & (width(bam_gr) <= 23)]

### count reads 
# overlap
overlap_reads_all <- GenomicRanges::countOverlaps(query = exons_gr, subject = bam_gr_filt, type = "any", ignore.strand = T)
overlap_reads_sense <- GenomicRanges::countOverlaps(query = exons_gr, subject = bam_gr_filt, type = "any", ignore.strand = F)
overlap_reads_antisense <- overlap_reads_all - overlap_reads_sense

# save count table
overlap_counts <- 
  tibble(gene_id = names(overlap_reads_all), 
         count_sense = overlap_reads_sense,
         count_antisense = overlap_reads_antisense, 
         sample = bam_name, 
         experiment = experiment_name) %>% 
  readr::write_csv(x = ., path = file.path(outpath, str_c("smallRNA.perfect.21to23.counts", experiment_name, bam_name, "csv", sep = ".")))

