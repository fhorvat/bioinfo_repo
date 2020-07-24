### INFO: creates table with relations between ensembl and USCS seqnames using assembly report from NCBI
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1/gtf_with_some_refSeq_genes")

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
# transform last column in .gtf from refSeq to Ensembl style
refSeqToEnsemblStyleGtf <- function(gene_name){
  
  # find gene in RefSeq gtf
  refseq_gtf_filt <- refseq_gtf[str_detect(mcols(refseq_gtf)$gene, gene_name) %>% replace(., is.na(.), F)]
  
  # create ensembl style last column
  if(length(refseq_gtf_filt) > 0){
    
    refseq_last_column <- DataFrame(source = "refSeq", 
                                    type = mcols(refseq_gtf_filt)$type, 
                                    score = mcols(refseq_gtf_filt)$score, 
                                    phase = mcols(refseq_gtf_filt)$phase, 
                                    gene_id = mcols(refseq_gtf_filt)$gene, 
                                    gene_version = 1,
                                    gene_name = mcols(refseq_gtf_filt)$gene, 
                                    gene_source = "refSeq", 
                                    gene_biotype = mcols(refseq_gtf_filt)$gene_biotype, 
                                    transcript_id = mcols(refseq_gtf_filt)$transcript_id, 
                                    transcript_version = ifelse(mcols(refseq_gtf_filt)$type != "gene", 1, NA), 
                                    transcript_name = mcols(refseq_gtf_filt)$transcript_id, 
                                    transcript_source = ifelse(mcols(refseq_gtf_filt)$type != "gene", "refSeq", NA), 
                                    transcript_biotype = ifelse(mcols(refseq_gtf_filt)$type != "gene", mcols(refseq_gtf_filt)$gene_biotype, NA), 
                                    exon_number = NA, 
                                    exon_id = NA, 
                                    exon_version = ifelse(mcols(refseq_gtf_filt)$type == "exon", 1, NA), 
                                    protein_id = mcols(refseq_gtf_filt)$protein_id, 
                                    protein_version = ifelse(mcols(refseq_gtf_filt)$type == "CDS", 1, NA))
    
    # change last column to ensembl style
    mcols(refseq_gtf_filt) <- refseq_last_column
    
    # change seqnames
    new_seqname <- scaffold_names[scaffold_names$refSeq_accn == as.character(seqnames(refseq_gtf_filt)[1]), ]$ensembl_name
    seqlevels(refseq_gtf_filt) <- unique(c(seqlevels(refseq_gtf_filt), seqlevels(ensembl_gtf)))
    seqnames(refseq_gtf_filt) <- factor(new_seqname, levels = seqlevels(refseq_gtf_filt))
    
    # return 
    return(refseq_gtf_filt)
    
  }else{
    
    # error
    stop("Gene name not present in refSeq .gtf")
    
  }
  
}

# get gene info from refSeq
refSeqToGeneInfo <- function(gene_name){
  
  # find gene in RefSeq gtf
  refseq_gtf_filt <- refseq_gtf[str_detect(mcols(refseq_gtf)$gene, gene_name) %>% replace(., is.na(.), F)]
  
  # create geneInfo table
  if(length(refseq_gtf_filt) > 0){
    
    # get gene
    refseq_gtf_gene <- refseq_gtf_filt[mcols(refseq_gtf_filt)$type == "gene"]
    
    # create tibble
    refseq_geneInfo <- tibble(gene_id = mcols(refseq_gtf_gene)$gene, 
                              seqnames = scaffold_names[scaffold_names$refSeq_accn == as.character(seqnames(refseq_gtf_gene)), ]$ensembl_name,
                              start = start(refseq_gtf_gene), 
                              end = end(refseq_gtf_gene), 
                              strand = as.character(strand(refseq_gtf_gene)), 
                              gene_name =  mcols(refseq_gtf_gene)$gene, 
                              gene_biotype = mcols(refseq_gtf_gene)$gene_biotype, 
                              gene_description = str_c(mcols(refseq_gtf_gene)$Dbxref, ", added from RefSeq .gtf"))
    
    
    # return 
    return(refseq_geneInfo)
    
  }else{
    
    # error
    stop("Gene name not present in refSeq .gtf")
    
  }
  
}

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# ensembl gtf path
ensembl_gtf_path <- file.path(inpath, "..", "ensembl.99.MesAur1.0.20200415.UCSCseqnames.gtf.gz")

# refSeq gtf path
refseq_gtf_path <- file.path(inpath, "../../vfranke", "GCF_000349665.1_MesAur1.0_genomic.gtf")

# gene info path
gene_info_path <- file.path(inpath, "..", "ensembl.99.MesAur1.0.20200415.UCSCseqnames.geneInfo.csv")

# scaffold names path
scaffold_names_path <- file.path(inpath, "..", "GCF_000349665.1_MesAur1.0.ensembl2UCSC.txt")

######################################################## READ DATA
# read gtf using rtracklayer
ensembl_gtf <- rtracklayer::import.gff(ensembl_gtf_path)
refseq_gtf <- rtracklayer::import.gff(refseq_gtf_path)

# read gene info
gene_info <- readr::read_csv(gene_info_path)

# read scaffold names
scaffold_names <- readr::read_delim(scaffold_names_path, delim = "\t")

######################################################## MAIN CODE
# create list of genes you want to add to Ensembl .gtf
gene_list <- c("Piwil3", "Tex101")

# create new gtf name
gene_list_collapsed <-str_c(gene_list, collapse = "_")
gtf_name <- ensembl_gtf_path %>% basename %>% str_replace("gtf\\.gz", str_c(gene_list_collapsed, ".from_RefSeq.gtf"))


### add gene(s) from refSeq to Ensembl
# add refSeq seqnames to Ensembl .gtf
seqlevels(ensembl_gtf) <- unique(c(seqlevels(refseq_gtf), seqlevels(ensembl_gtf)))

# get genes from refSeq
refseq_genes_gtf <- 
  purrr::map(gene_list, refSeqToEnsemblStyleGtf) %>% 
  do.call(c, .)

# add genes 
ensembl_gtf <- c(ensembl_gtf, refseq_genes_gtf)

# save
rtracklayer::export.gff(ensembl_gtf, file.path(outpath, gtf_name))

# gzip
system(stringr::str_c("gzip ", file.path(outpath, gtf_name)))


### add info about genes to gene info table
# create table
refseq_genes_info <- 
  purrr::map(gene_list, refSeqToGeneInfo) %>% 
  dplyr::bind_rows(.)

# bind with gene info
gene_info <- bind_rows(gene_info, refseq_genes_info)

# save
readr::write_csv(gene_info, file.path(outpath, str_replace(gtf_name, "gtf", "geneInfo.csv")))
