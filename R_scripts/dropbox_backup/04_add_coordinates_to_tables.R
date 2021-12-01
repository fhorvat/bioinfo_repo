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

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
# library(AnnotationDbi)
# library(goseq)
# 
# library(org.Mm.eg.db)
# library(EnsDb.Mmusculus.v79)
# library(GO.db)

######################################################## PATH VARIABLES
outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/diffExp/results"
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/diffExp/results"
gtf_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Ensembl_GRCm38.86.20161128.gtf.gz"
diff_list <- list.files(inpath, pattern = "diffExp_significant.csv", full.names = T)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "GffToGRanges.R"))
source(file.path(lib_path, "headt.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
ensembl_gtf <- read_delim(file = gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c")) 

######################################################## MAIN CODE
# get exons for each transcript
gtf_trans <- GffToGRanges(ensembl_gtf, "exon")

# get unique gene_id/transcript_id combinations
gids <- unique(values(gtf_trans)[c('gene_id', 'transcript_id')])

# splits exon ranges on transcripts
gtf_trans <- split(gtf_trans, gtf_trans$transcript_id)

# orders transcripts based on number of exons in transcripts
gtf_trans <- gtf_trans[order(elementNROWS(gtf_trans), decreasing = T)] 

# keeps only first transcript of each gene (the one with most exons)
gtf_trans <- gtf_trans[!duplicated(gids$gene_id[match(names(gtf_trans), gids$transcript_id)])]

# get the full range of transcripts
gtf_trans_full <-
  range(gtf_trans) %>%  
  unlist() %>% 
  as.character() %>%  
  tibble(transcript_id = names(.), full_pos = .)
  
# read differential expressed genes tables
diff_table_list <- 
  lapply(diff_list, readr::read_csv) %>% 
  magrittr::set_names(basename(diff_list) %>% 
                        stringr::str_replace(., ".csv", ""))

# join coordinates on all list
diff_table_list_clean <- 
  lapply(X = names(diff_table_list), FUN = function(X){
    
    # join coordinates
    left_join(diff_table_list[[X]], gtf_trans_full, by = "transcript_id") %>% 
      dplyr::select(transcript_id, full_pos, everything()) %T>% 
      readr::write_csv(., path = file.path(outpath, str_c(X, "_with_position.csv")))
  
})


