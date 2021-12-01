### INFO: 
### DATE: 08. 09. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Data/Mapped/STARlong_mm10/unspliced_reads")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tibble)
library(data.table)

library(Rsamtools)
library(Biostrings)
library(GenomicRanges)
library(GenomicAlignments)

######################################################## PATH VARIABLES
outpath <- getwd()

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "bamToGRangesList.R"))
source(file.path(lib_path, "vfranke", "GffToGRanges.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
# ENSEMBL gtf, convert to GRanges
gtf_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.gtf.gz"
gtf <- read.table(gtf_path, header = FALSE, sep = "\t", stringsAsFactors = F)
gtf[, 1] <- paste0("chr", gtf[, 1])
gtf <- gtf[!str_detect(gtf[, 1], "NT"),]
gtf[gtf[, 1] == "chrMT", 1] <- "chrM"
gtf <- GffToGRanges(gtf, "exon")

######################################################## MAIN CODE
##### for each gene get all exons
gtf_genes <-
  as.data.frame(gtf) %>%
  as.tibble() %>%
  dplyr::distinct(exon_id, .keep_all = T) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
  split(., .$gene_id)

##### for each gene get transcript with most exons
# create ENSEMBL transcripts ranges, takes one transcript for each gene (the one with the most exons)
gids <- unique(values(gtf)[c("gene_id", "transcript_id")])
gtf_trans <- split(gtf, gtf$transcript_id)
gtf_trans <- gtf_trans[order(-elementNROWS(gtf_trans))]
gtf_trans <- gtf_trans[!duplicated(gids$gene_id[match(names(gtf_trans), gids$transcript_id)])]

# take transcripts with single exon
gtf_trans_single <- gtf_trans[elementNROWS(gtf_trans) == 1]
gtf_trans_single <- unlist(gtf_trans_single)
gtf_trans_single$ex.num <- 1
gtf_trans_single$ex.tot <- 1
gtf_trans_single <- split(gtf_trans_single, gtf_trans_single$transcript_id)
gtf_trans_single <- gtf_trans_single[elementNROWS(gtf_trans_single) != 0]

# take transcripts with more than one exon, count and enumerate exons in each transcript
gtf_trans <- gtf_trans[elementNROWS(gtf_trans) > 1]
gtf_trans <- unlist(gtf_trans)
d_val <- data.table(as.data.frame(values(gtf_trans)))
d_val$strand <- as.character(strand(gtf_trans))
d_val[d_val$strand == "+" , `:=`( COUNT = .N , IDX = 1:.N ) , by = transcript_id[strand == "+"]]
d_val[d_val$strand == "-" , `:=`( COUNT = .N , IDX = .N:.1) , by = transcript_id[strand == "-"]]
gtf_trans$ex.num <- d_val$IDX
gtf_trans$ex.tot <- d_val$COUNT
gtf_trans <- split(gtf_trans, as.character(gtf_trans$transcript_id))
gtf_trans <- c(gtf_trans, gtf_trans_single)

##### for each transcript get all exons
gtf_exons <- split(gtf, as.character(gtf$transcript_id))

##### for each transcript get all introns
gtf_introns <- 
  gtf_exons %>%
  as.data.frame(., row.names = NULL) %>% 
  as.tibble(.) %>% 
  dplyr::select(seqnames, start, end, strand, transcript_id, start_exon = exon_id) %>%
  dplyr::group_by(transcript_id) %>% 
  dplyr::arrange(start) %>% 
  dplyr::mutate(start_int = end + 1,
                end_int = lead(start) - 1) %>%  
  dplyr::filter(!is.na(end_int)) %>% 
  dplyr::mutate(start = ifelse(strand == "-", -(start), start)) %>%
  dplyr::arrange(start) %>%
  dplyr::mutate(start = abs(start)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-start, -end) %>% 
  dplyr::select(seqnames, start = start_int, end = end_int, everything()) %>% 
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T) %>% 
  split(., as.character(.$transcript_id))
