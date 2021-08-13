### INFO: 
### DATE: Sat Aug 17 18:47:14 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/arrays/Joshi_2007_BMCDevBiol_GSE5558")

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
library(GEOquery)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- file.path(inpath, "align_probes/minimap2")

# tidy annotation table path
annotation_tidy_path <- file.path(inpath, "GPL3771.annot.tidy.csv")

# gene info path
gene_info_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.geneInfo.csv"

# mapped probes path
mapped_probes_path <- file.path(inpath, "align_probes", "GPL3771.annot.GenBank.mm10.bam")

# probe fasta index path
fai_path <- file.path(inpath, "align_probes", "GPL3771.annot.GenBank.fasta.fai")
  
######################################################## READ DATA
# read annotation table
annotation_tidy <- readr::read_csv(annotation_tidy_path)

# read basic annotation 
annotation_all <- getGEO("GPL3771", destdir = inpath, AnnotGPL = T)

# read gene info
gene_info <- readr::read_csv(gene_info_path)

# read probes fasta index
fai <- readr::read_delim(fai_path, delim = "\t", col_names = F)

# read mapped probes coordinates
mapped_probes <- GenomicAlignments::readGAlignmentsList(mapped_probes_path, use.names = T, 
                                                        param = ScanBamParam(what = "qname", tag = "NM"))

######################################################## MAIN CODE
### which probes are not on GenBank?
annotation_absent <- 
  annotation_all %>% 
  Table(.) %>% 
  as_tibble(.) %>% 
  dplyr::filter(!(`GenBank Accession` %in% fai$X1)) %>% 
  readr::write_csv(., path = file.path(outpath, "probes_not_in_GenBank.csv"))

### prepare data
# get genes coordinates
gene_gr <- 
  gene_info %>% 
  GenomicRanges::GRanges(.)

### tidy probe alignments
# get probe alignments
probe_align <- 
  mapped_probes %>% 
  unlist(.)

# set unique names to alignments
names(probe_align) <- make.unique(names(probe_align))

# get ranges of only aligned part
probe_ranges <- 
  probe_align %>% 
  GenomicRanges::grglist(.)


### find which probe alignes to which gene
# find overlaps between probes and alignments
probe_gene_overlaps <- findOverlaps(probe_ranges, gene_gr, ignore.strand = T)

# extract probes and gene info
probe_gene_tb <- 
  tibble(probe = probe_ranges[queryHits(probe_gene_overlaps)] %>% names, 
                        gene_name = gene_gr[subjectHits(probe_gene_overlaps)]$gene_name) %>% 
  dplyr::mutate(probe = str_remove(probe, "\\.[0-9]+")) %>% 
  dplyr::group_by(probe) %>% 
  dplyr::summarise(gene_name = str_c(gene_name, collapse = ".")) %>% 
  dplyr::ungroup(.) %>% 
  cbind(., setDT(.)[, tstrsplit(gene_name, "\\.")]) %>% 
  as_tibble(.) %>% 
  dplyr::select(GenBank = probe, gene_name, V2, V3)

# join with annotation table
annotation_tidy_aligned <- 
  annotation_tidy %>% 
  dplyr::rename(gene_description.original = gene_name) %>% 
  left_join(., probe_gene_tb, by = "GenBank")



