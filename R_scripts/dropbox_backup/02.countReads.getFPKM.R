#!/common/WORK/fhorvat/programi/R/R-3.4.3/bin/Rscript
### RUN: qsub -q MASTER -M fihorvat@gmail.com -m n -N pbs.Rscript.summarizeOverlaps -l select=ncpus=10:mem=300g -j oe 01.expression_analysis.R
### INFO: Count reads and get FPKM of genes in mESC and oocytes sequenced in February 2018
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat

rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Analysis/expression")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "GffToGRanges.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set output path
outpath <- getwd()

# sample table
sample_path <- "/common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Data/Documentation/mESC_oocytes_2018.sample_table.csv"

# reduced exons path
exons_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.91.GRCm38.p5.20180512.UCSCseqnames.reducedExons.RDS"

# info about genes path
genes_info_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.91.GRCm38.p5.20180512.UCSCseqnames.geneInfo.csv"

######################################################## READ DATA
# read sample table
sample_table <- readr::read_csv(file = sample_path)

# read ENSEMBL reduced exons 
exons_gr <- readRDS(file = exons_path)

# read additional info about genes
genes_info <- readr::read_csv(genes_info_path)

######################################################## MAIN CODE
# get total length of all exons for each transcript
exons_width <- 
  width(exons_gr) %>%
  sum(.) %>% 
  tibble(gene_id = names(.), width = .)

# # get count of reads, save summarizedExperiment as RDS
# bamfiles <- Rsamtools::BamFileList(sample_table$bam_path, yieldSize = 2000000)
# BiocParallel::register(BiocParallel::MulticoreParam())
# se <- GenomicAlignments::summarizeOverlaps(features = exons_gr,
#                                            reads = bamfiles,
#                                            mode = "Union",
#                                            singleEnd = TRUE,
#                                            ignore.strand = TRUE)
# 
# saveRDS(se, file = file.path(outpath, "GRCm38.91.reducedExons.summarizedOverlaps.RDS"))

# read summarizedExperiment from RDS file
se <- readRDS(file = file.path(outpath, "GRCm38.91.reducedExons.summarizedOverlaps.RDS")) 

### FPKM
# get data.frame of counts, transform to FPKM
fpkm_df <-
  assay(se) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  as.tibble(.) %>%
  set_colnames(., str_replace(colnames(.), ".genome.Aligned.sortedByCoord.out.bam", "")) %>%
  tidyr::gather(key = sample_id, value = counts, -gene_id) %>%
  dplyr::left_join(., sample_table %>% dplyr::select(sample_id, library_size), by = "sample_id") %>%
  dplyr::left_join(., exons_width, by = "gene_id") %>%
  dplyr::mutate(library_size = round(library_size / 1E6, 6),
                width = round(width / 1E3, 3),
                fpm = (counts / library_size),
                fpkm = (fpm / width)) %>%
  dplyr::select(gene_id, sample_id, fpkm) %>%
  tidyr::spread(key = sample_id, value = fpkm) %>%
  dplyr::left_join(., genes_info, by = "gene_id") %T>%
  readr::write_csv(., path = file.path(outpath, "GRCm38.91.reducedExons.FPKM.csv"))
