### INFO: creates table with relations between ensembl and USCS seqnames using assembly report from NCBI
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1/RefSeq")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicRanges)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# refSeq gtf path
refseq_gtf_path <- file.path(inpath, "GCF_000349665.1_MesAur1.0_genomic.gtf.gz")

# scaffold names path
scaffold_names_path <- file.path(inpath, "GCF_000349665.1_MesAur1.0.ensembl2UCSC.txt")

######################################################## READ DATA
# read gtf using rtracklayer
refseq_gtf <- rtracklayer::import.gff(refseq_gtf_path)

# read scaffold names
scaffold_names <- readr::read_delim(scaffold_names_path, delim = "\t")

######################################################## MAIN CODE
# set GTF name
gtf_name <- refseq_gtf_path %>% basename %>% str_replace("gtf\\.gz", "UCSCseqnames.gtf")

# add ENSEMBL scaffold names
refseq_gtf_ensembl_seqnames <- 
  refseq_gtf %>% 
  as_tibble(.) %>% 
  dplyr::left_join(., scaffold_names %>% dplyr::select(seqnames = refSeq_accn, ensembl_seqnames = ensembl_name), by = "seqnames") %>% 
  dplyr::select(-seqnames) %>% 
  dplyr::select(seqnames = ensembl_seqnames, everything()) %>% 
  dplyr::mutate(product = replace(product, gene_id == "", make.unique(product[gene_id == ""])), 
                gene_id = ifelse(gene_id == "", product, gene_id)) %>% 
  GRanges(.)

# save
rtracklayer::export.gff(refseq_gtf_ensembl_seqnames, file.path(outpath, gtf_name))

# gzip
system(stringr::str_c("gzip ", file.path(outpath, gtf_name)))


### get gene info
# get gene
refseq_gtf_gene <- refseq_gtf_ensembl_seqnames[mcols(refseq_gtf_ensembl_seqnames)$type == "gene"]

# create tibble
refseq_geneInfo <- tibble(gene_id = mcols(refseq_gtf_gene)$gene, 
                          seqnames = as.character(seqnames(refseq_gtf_gene)),
                          start = start(refseq_gtf_gene), 
                          end = end(refseq_gtf_gene), 
                          strand = as.character(strand(refseq_gtf_gene)), 
                          gene_name = mcols(refseq_gtf_gene)$gene, 
                          gene_biotype = mcols(refseq_gtf_gene)$gene_biotype, 
                          gene_description = mcols(refseq_gtf_gene)$db_xref)

# save
readr::write_csv(refseq_geneInfo, file.path(outpath, str_replace(gtf_name, "gtf", "geneInfo.csv")))




