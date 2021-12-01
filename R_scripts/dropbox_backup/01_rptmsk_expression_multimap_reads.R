### INFO: expression of multimap reads (> 10 mapping location)
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

# library(org.Mm.eg.db)
# library(EnsDb.Mmusculus.v79)
# library(GO.db)
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)

######################################################## PATH VARIABLES
outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/STAR_mm10_noMultimapFilter/multimap_reads/analysis"
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/STAR_mm10_noMultimapFilter/multimap_reads"
bam_list <- list.files(path = bam_path, pattern = "*bam$", recursive = F, full.names = T)
log_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/STAR_mm10_noMultimapFilter"
log_list <- list.files(path = log_path, pattern = "*Log.final.out$", recursive = T, full.names = T)

lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
repeatmasker_path <- "/common/WORK/fhorvat/reference/mouse/mm10/UCSC/UCSC_repeatMaskerVIZ_OutBaseline_20160824.txt.gz"
ensembl_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.gtf.gz"

######################################################## SOURCE FILES
source(file.path(lib_path, "GffToGRanges.R"))
source(file.path(lib_path, "headt.R"))
source("/common/WORK/fhorvat/code_library/R_scripts/vfranke/BamWorkers.R")

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

# repeatMaskerVIZ 
rptmsk_gr <-
  read_delim(file = repeatmasker_path, delim = "\t") %>%
  dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, repName, repClass, repFamily) %>%
  dplyr::mutate(full_pos = str_c(seqnames, ":", start, "-", end, "|", strand, "|", repClass, "|", repName)) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

# get length of rptmsk
rptmsk_width <-
  tibble(full_pos = rptmsk_gr$full_pos, 
         width = width(rptmsk_gr))

######################################################## MAIN CODE
# ### take transcript with most exons for each gene
# # ENSEMBL gtf
# ensembl_gtf <- read_delim(file = ensembl_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c")) 
# gtf_trans <- GffToGRanges(ensembl_gtf, "exon")
# 
# # get total length of all exons for each transcript
# trans_width <- 
#   sum(width(gtf_trans)) %>% 
#   tibble(transcript_id = names(.), width = .)
# 
# # get GRanges from .gtf, filter only exons
# gtf_trans <- gtf_trans[gtf_trans$gene_biotype %in% c("protein_coding", "lincRNA")]
# 
# # get unique gene_id/transcript_id combinations
# gids <- unique(values(gtf_trans)[c('gene_id', 'transcript_id')])
# 
# # splits exon ranges on transcripts
# gtf_trans <- split(gtf_trans, gtf_trans$transcript_id)
# 
# # orders transcripts based on number of exons in transcripts
# gtf_trans <- gtf_trans[order(elementNROWS(gtf_trans), decreasing = T)] 
# 
# # keeps only first transcript of each gene (the one with most exons)
# gtf_trans <- gtf_trans[!duplicated(gids$gene_id[match(names(gtf_trans), gids$transcript_id)])]
# 
# 
# ### get count of reads over transcripts
# # ensembl 
# bamfiles <- Rsamtools::BamFileList(bam_list, yieldSize = 2000000)
# BiocParallel::register(BiocParallel::MulticoreParam())
# se_ensembl <- GenomicAlignments::summarizeOverlaps(features = gtf_trans,
#                                                    reads = bamfiles,
#                                                    mode = "Union",
#                                                    singleEnd = FALSE,
#                                                    ignore.strand = TRUE)
# saveRDS(se_ensembl, file = file.path(outpath, "se_ensembl_multimap_reads.RDS"))
# 
# # ensembl FPKM df
# fpkm_df <- 
#   assay(se_ensembl) %>% 
#   as.tibble(.) %>% 
#   magrittr::set_colnames(., value = stringr::str_replace(colnames(.), "_multimap.bam", "")) %>% 
#   dplyr::mutate(transcript_id = names(gtf_trans)) %>% 
#   dplyr::mutate(count_sum = rowSums(.[-ncol(.)])) %>% 
#   dplyr::filter(count_sum > 0) %>% 
#   dplyr::select(-count_sum) %>% 
#   tidyr::gather(key = sample, value = counts, -transcript_id) %>% 
#   dplyr::left_join(., sample_table, by = "sample") %>% 
#   dplyr::left_join(., trans_width, by = "transcript_id") %>% 
#   dplyr::mutate(library_size = round(library_size / 1E6, 2), 
#                 width = round(width / 1E3, 2), 
#                 fpm = (counts / library_size), 
#                 fpkm = (fpm / width)) %>% 
#   dplyr::select(transcript_id, sample, fpkm) %>% 
#   tidyr::spread(key = sample, value = fpkm) %>% 
#   dplyr::mutate(gene_symbol = mapIds(EnsDb.Mmusculus.v79, keys = transcript_id, column = "SYMBOL", keytype = "TXID", multiVals = "first"),
#                 gene_biotype = mapIds(EnsDb.Mmusculus.v79, keys = transcript_id, column = "GENEBIOTYPE", keytype = "TXID", multiVals = "first"))
# 
# gene_name <- mapIds(org.Mm.eg.db, keys = fpkm_df$gene_symbol, column = "GENENAME", keytype = "SYMBOL", multiVals = "first")
# gene_name[sapply(gene_name, is.null)] <- NA
# gene_name <- unname(unlist(gene_name))
# fpkm_df$gene_name <- gene_name
# 
# readr::write_csv(x = fpkm_df, path = file.path(outpath, "FPKM_ensembl_multimap_reads.csv"))

# # repeatMasker 
# se_rptmsk <- GenomicAlignments::summarizeOverlaps(features = rptmsk_gr,
#                                                   reads = bamfiles,
#                                                   mode = "Union",
#                                                   singleEnd = FALSE,
#                                                   ignore.strand = TRUE)
# saveRDS(se_rptmsk, file = file.path(outpath, "se_rptmsk_multimap_reads.RDS"))
se_rptmsk <- readRDS(file = file.path(outpath, "se_rptmsk_multimap_reads.RDS"))
  
# repeatMasker FPKM df
fpkm_df <- 
  assay(se_rptmsk) %>% 
  as.tibble(.) %>% 
  magrittr::set_colnames(., value = stringr::str_replace(colnames(.), ".bam", "")) %>% 
  dplyr::mutate(full_pos = rptmsk_gr$full_pos) %>% 
  dplyr::mutate(count_sum = rowSums(.[-ncol(.)])) %>% 
  dplyr::filter(count_sum > 0) %>% 
  dplyr::select(-count_sum) %>% 
  tidyr::gather(key = sample, value = counts, -full_pos) %>% 
  dplyr::left_join(., sample_table, by = "sample") %>% 
  dplyr::left_join(., rptmsk_width, by = "full_pos") %>% 
  dplyr::mutate(library_size = round(library_size / 1E6, 2), 
                width = round(width / 1E3, 2), 
                fpm = (counts / library_size), 
                fpkm = round((fpm / width), 3)) %>% 
  dplyr::select(full_pos, sample, fpkm) %>% 
  tidyr::spread(key = sample, value = fpkm) %T>% 
  readr::write_csv(., path = file.path(outpath, "FPKM_rptmsk_multimap_reads.csv"))
  
  
