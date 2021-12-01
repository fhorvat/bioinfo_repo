### INFO: read smallRNA seq bam file, get counts over exons
### DATE: Mon May 14 12:21:35 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# mouse
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/mouse")
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# cow
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/cow")
genome_path <- "/common/DB/genome_reference/cow/bosTau8.UMD3.1.GCA_000003055.4"

# pig
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/pig")
genome_path <- "/common/DB/genome_reference/pig/susScr11.Sscrofa11.1.GCA_000003025.6"

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
# setdiff between two GRanges objects in a way which keeps names of final overlap
GRanges.parallelSetDiff <- function(gr1, gr2, ignore.strand = T){
  
  # find overlaps between two GRanges
  hits <- findOverlaps(gr1, gr2, ignore.strand = ignore.strand)
  
  # extract all overalaping features from subject as list
  grl <- extractList(gr2, as(hits, "List"))
  
  # parallel set difference - query vs. subject
  diff_list <- psetdiff(gr1, grl, ignore.strand = ignore.strand)
  
  # set names to diff. list, unlist
  names(diff_list) <- names(gr1)
  diff_gr <- unlist(diff_list)
  
}

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# get gene info path
info_path <- list.files(genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.geneInfo.csv", full.names = T)

# get paths of reduced exons
exons_path <- list.files(genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.reducedExons.RDS", full.names = T)

# get repeatMasker path
rmsk_path <- list.files(genome_path, pattern = "rmsk.*[0-9]{6}.clean.fa.out.gz", full.names = T)

# # get path to miRBase
# mir_path <- file.path(genome_path, "miRBase.22.mm10.20181605.gff3")

######################################################## READ DATA
# read info about genes
genes_info <- readr::read_csv(info_path)

# gtf exons by genes
exons_gr <- readRDS(file = exons_path)

# read repeatMasker
rmsk_df <- readr::read_delim(rmsk_path, delim = "\t")

# # read miRBase
# mir_df <- ape::read.gff(file = mir_path) 

######################################################## MAIN CODE
### prepare ranges to count over - exons with repeats and miRNA removed
# get GRanges from repeatMasker
rmsk_gr <- 
  rmsk_df %>% 
  GenomicRanges::GRanges(.) %>% 
  GenomicRanges::reduce(.)

# # get coordinates of miRNAs
# mir_gr <- 
#   mir_df %>% 
#   dplyr::select(seqnames = seqid, start, end, strand, type, attributes) %>% 
#   GenomicRanges::GRanges(.)

# micro RNA from ensembl
mir_gene_id <- 
  genes_info %>% 
  dplyr::filter(gene_biotype == "miRNA") %$% 
  gene_id

# #####
# # micro RNA from ensembl
# rrna <- 
#   genes_info %>% 
#   dplyr::filter(gene_biotype == "rRNA") %$% 
#   gene_id
# exons_rrna <- exons_gr[(names(exons_gr) %in% rrna)]
# 
# rmsk_rrna <- subsetByOverlaps(exons_rrna, rmsk_gr, invert = T)
# #####

### filter exons
# remove microRNA from ensembl
exons_gr <- exons_gr[!(names(exons_gr) %in% mir_gene_id)]

# substract ranges overlaping rmsk table
exons_gr <- 
  exons_gr %>% 
  unlist(.) %>% 
  GRanges.parallelSetDiff(gr1 = ., gr2 = rmsk_gr, ignore.strand = T) %>% 
  # GRanges.parallelSetDiff(gr1 = ., gr2 = mir_gr, ignore.strand = T) %>% 
  GenomicRanges::split(., names(.))

# save filtered exon ranges
saveRDS(object = exons_gr, file = file.path(outpath, "ensembl.91.Sscrofa11.1.20180404.UCSCseqnames.reducedExons.filt.rmsk_miRNA.RDS"))

