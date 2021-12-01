### INFO: 
### DATE: Wed Sep 09 09:23:26 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/otherSeq/Gahurova_2017_EpigeneticsChromatin_GSE86297/Analysis/methylation")

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

library(methylKit)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### genomic files
# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)


### experiment path
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/otherSeq/Gahurova_2017_EpigeneticsChromatin_GSE86297"

# sample table path
sample_tb_path <- file.path(base_path, "Data/Documentation")
sample_tb_path <- list.files(sample_tb_path, ".*sampleTable\\.csv$", full.names = T)

# list bam files
cpg_path_list <- file.path(base_path, "Data/Mapped/Bismark_mm10")
cpg_path_list <- list.files(cpg_path_list, ".*\\.bam$", full.names = T)[1:4]

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read sample table
sample_tb <- readr::read_csv(sample_tb_path)

######################################################## MAIN CODE
# read CpG methylation calls from Bismark bam
cpg_calls <- processBismarkAln(location = as.list(cpg_path_list),
                               sample.id = as.list(cpg_path_list %>% basename(.) %>% str_remove(., "_bismark_bt2\\.bam")), 
                               treatment = c(1, 1, 0, 0), 
                               assembly = "mm10", 
                               read.context = "CpG", 
                               save.folder = outpath)

# get methylation stats
pdf(file = file.path(outpath, str_c("expl_plot", "methylation_distribution", "pdf", sep = ".")))
for(i in 1:length(cpg_calls)){getMethylationStats(cpg_calls[[i]], plot = TRUE, both.strands = FALSE)}
dev.off()

# get per base information
pdf(file = file.path(outpath, str_c("expl_plot", "read_coverage_per_base", "pdf", sep = ".")))
for(i in 1:length(cpg_calls)){getCoverageStats(cpg_calls[[i]], plot = TRUE, both.strands = FALSE)}
dev.off()

# unite samples
cpg_calls_united <- unite(cpg_calls, destrand = FALSE)

# get sample correlation
png(file = file.path(outpath, str_c("expl_plot", "correlation", "png", sep = ".")))
getCorrelation(cpg_calls_united, plot = TRUE)
dev.off()
