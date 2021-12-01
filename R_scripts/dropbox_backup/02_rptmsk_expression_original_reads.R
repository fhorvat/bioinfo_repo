#!/home/students/fhorvat/R/bin/Rscript
### qsub -V -q MASTER -l select=8:mem=200gb -N sumOverlaps -j oe -o /common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/STAR_mm10_noMultimapFilter/multimap_reads/analysis rptmsk_expression_original_reads_RDS.R
### INFO: expression of original reads (<= 10 mapping location)
### DATE: 25. 08. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()
# options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/STAR_mm10_noMultimapFilter/multimap_reads/analysis")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tibble)
library(data.table)
library(doMC)

library(Rsamtools)
library(Biostrings)
library(GenomicRanges)
library(GenomicAlignments)
library(BiocParallel)

######################################################## PATH VARIABLES
outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/STAR_mm10_noMultimapFilter/multimap_reads/analysis"
bam_path <- "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2"
bam_list <- list.files(path = bam_path, pattern = "*bam$", recursive = T, full.names = T)
log_list <- list.files(path = bam_path, pattern = "*Log.final.out$", recursive = T, full.names = T)
multimap_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/STAR_mm10_noMultimapFilter/multimap_reads/analysis/FPKM_rptmsk_multimap_reads.csv"

lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
repeatmasker_path <- "/common/WORK/fhorvat/reference/mouse/mm10/UCSC/UCSC_repeatMaskerVIZ_OutBaseline_20160824.txt.gz"
ensembl_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.gtf.gz"

######################################################## SOURCE FILES
source(file.path(lib_path, "GffToGRanges.R"))
source(file.path(lib_path, "headt.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
# experiment table
sample_table <- 
  tibble(sample = str_replace_all(bam_list, "\\/.*\\/|.bam", ""), 
         log_path = log_list) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(library_size = readr::read_lines(file = log_path) %>% 
                  str_subset(., "Uniquely mapped reads number|Number of reads mapped to multiple loci") %>% 
                  str_extract(., "[0-9].*") %>% 
                  as.integer(.) %>% 
                  sum(.)) %>% 
  dplyr::select(-log_path)

# multimap reads features
rptmsk_feature <- 
  read_csv(file = multimap_path) %$%
  full_pos
  
# repeatMaskerVIZ - retrotransposons only
rptmsk_gr <-
  read_delim(file = repeatmasker_path, delim = "\t") %>%
  dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, repName, repClass, repFamily) %>%
  dplyr::mutate(full_pos = str_c(seqnames, ":", start, "-", end, "|", strand, "|", repClass, "|", repName)) %>%
  dplyr::filter(full_pos %in% rptmsk_feature) %>% 
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

# get length of rptmsk
rptmsk_width <-
  tibble(full_pos = rptmsk_gr$full_pos, 
         width = width(rptmsk_gr))

######################################################## MAIN CODE
# ### get count of reads over transcripts
# # repeatMasker
bamfiles <- Rsamtools::BamFileList(bam_list, yieldSize = 2000000)
BiocParallel::register(BiocParallel::MulticoreParam())
se_rptmsk <- GenomicAlignments::summarizeOverlaps(features = rptmsk_gr,
                                                  reads = bamfiles,
                                                  mode = "Union",
                                                  singleEnd = FALSE,
                                                  ignore.strand = TRUE)
# saveRDS(se_rptmsk, file = file.path(outpath, "se_rptmsk_original_reads.RDS"))
se_rptmsk <- readRDS(file = file.path(outpath, "se_rptmsk_original_reads.RDS"))

# repeatMasker FPKM df
fpkm_df <-
  assay(se_rptmsk) %>%
  as.tibble(.) %>%
  magrittr::set_colnames(., value = stringr::str_replace(colnames(.), ".bam", "")) %>%
  dplyr::mutate(full_pos = rptmsk_gr$full_pos) %>%
  tidyr::gather(key = sample, value = counts, -full_pos) %>%
  dplyr::left_join(., sample_table, by = "sample") %>%
  dplyr::left_join(., rptmsk_width, by = "full_pos") %>%
  dplyr::mutate(library_size = round(library_size / 1E6, 2),
                width = round(width / 1E3, 2),
                fpm = (counts / library_size),
                fpkm = round((fpm / width), 3)) %>%
  dplyr::select(full_pos, sample, fpkm) %>%
  tidyr::spread(key = sample, value = fpkm) %T>%
  readr::write_csv(., path = file.path(outpath, "FPKM_rptmsk_original_reads.csv"))


