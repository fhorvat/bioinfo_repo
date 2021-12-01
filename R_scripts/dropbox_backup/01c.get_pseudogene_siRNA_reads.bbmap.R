### INFO: 
### DATE: Sun Nov 24 12:22:19 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/siRNA_pseudogenes")

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

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(DESeq2)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)

# bam paths
bam_path <- list.files(inpath, ".*bam")

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read bam
bam_gr <- GenomicAlignments::readGAlignments(file = bam_path[2], param = ScanBamParam(tag = c("nM", "NH"), what = "qname"))

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

######################################################## MAIN CODE
# get coordinates of Sirena1 pseudogene (Gm31941), Elob and Elobl
genes_gr <- 
  genes_info %>% 
  dplyr::filter(gene_name %in% c("Gm31941", "Elob", "Elobl")) %>% 
  dplyr::select(gene_id, gene_name, seqnames:strand) %>% 
  # dplyr::mutate(start = replace(start, gene_name == "Gm31941",  46975357),
  #               end = replace(end, gene_name == "Gm31941",  46975758)) %>%
  GRanges(.)

# save as bed
rtracklayer::export.bed(genes_gr, file.path(outpath, "Sirena1_pseudo.Elob.Elobl.bed"))

### get sequences and write to separate fasta files
# get genomic sequences
genes_fasta <- Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm10, genes_gr)
names(genes_fasta) <- mcols(genes_gr)$gene_name

# write as one fasta
Biostrings::writeXStringSet(x = genes_fasta, filepath = file.path(outpath, "Sirena1_pseudo.Elob.Elobl.fa"))

# write in loop
for(fasta in names(genes_fasta)){
  
  # get one sequence and write
  genes_fasta[fasta] %>%
    writeXStringSet(., filepath = file.path(outpath, str_c("gene_seq", fasta, "fa", sep = ".")), format = "fasta")
  
}


### get reads mapping to our genes of interest
# find overlaps
overlaps <- findOverlaps(query = bam_gr, subject = genes_gr, ignore.strand = T, type = "any")

# subset bam
bam_hits <- 
  bam_gr[queryHits(overlaps)] %>% 
  as_tibble(.) %>% 
  mutate(gene_name = mcols(genes_gr[subjectHits(overlaps)])$gene_name) %>% 
  group_by(qname) %>% 
  dplyr::summarise(gene_name = str_c(gene_name, collapse = ","), 
                   nM = str_c(nM, collapse = ","), 
                   NH = str_c(NH, collapse = ","))



# get all alignments of those reads
bam_all <- bam_gr[which(mcols(bam_gr)$qname %in% bam_hits$qname)]

# annotate
overlaps <- findOverlaps(bam_all, exons_gr, ignore.strand = T)

# add info about genes
bam_annotate <- 
  bam_all[queryHits(overlaps)] %>% 
  as_tibble(.) %>% 
  dplyr::mutate(gene_id = names(exons_gr[subjectHits(overlaps)])) %>% 
  dplyr::left_join(., genes_info %>% dplyr::select(gene_id, gene_name), by = "gene_id") %>% 
  group_by(qname) %>% 
  dplyr::summarise(gene_name = str_c(gene_name, collapse = ", "), 
                   nM = str_c(nM, collapse = ", "), 
                   NH = str_c(unique(NH), collapse = ", "))

