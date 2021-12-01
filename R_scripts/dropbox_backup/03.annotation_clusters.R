### INFO:
### DATE: Tue Jan 07 15:36:14 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/mESC_MosIR/datasets/Jul_2018/Analysis/expression/smallRNA_clusters")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(data.table)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
### scans chromosomes in bam file
chrFinder <- function(bam.path, filter = FALSE, output = "data.frame") {
  
  # scan bam
  s <- Rsamtools::scanBamHeader(bam.path)
  st <- s[[1]]$text
  st <- do.call(rbind, st[names(st) == "@SQ"])
  st[, 1] <- str_replace(st[, 1], "SN:", "")
  st[, 2] <- str_replace(st[, 2], "LN:", "")
  
  # filter
  if (filter == TRUE) {
    st <- st[!str_detect(st[, 1], "random")]
  }
  
  # output
  if (output == 'data.frame') {
    vst <- data.frame(chr = st[, 1],
                      chrlen = as.numeric(st[, 2]),
                      stringsAsFactors = F)
    
  } else{
    vst <- as.numeric(st[, 2])
    names(vst) <- st[, 1]
    
  }
  
  # return
  return(vst)
  
}

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_path <-"/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# get gene info path
info_path <- file.path(genome_path, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.geneInfo.csv")

# get paths of reduced exons
exons_path <- file.path(genome_path, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.reducedExons.RDS")

# get repeatMasker path
rmsk_path <- file.path(genome_path, "rmsk.mm10.20180919.clean.fa.out.gz")

# get miRBase gff path
mirbase_path <- file.path(genome_path, "miRBase.22.mm10.20181605.gff3")

# get path to sample mapped bam file
sample_bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/mESC_MosIR/datasets/Jul_2018/Data/Mapped/STAR_mm10/2_original_mapping/s_RS10_NC_r1.SE.bam"

######################################################## READ DATA
# read info about genes
genes_info <- readr::read_csv(info_path)

# gtf exons by genes
exons_gr <- readRDS(file = exons_path)

# read repeatMasker
rmsk_df <- readr::read_delim(rmsk_path, delim = "\t")

# read miRBase gff
mirna_gr <- rtracklayer::import.gff(con = mirbase_path)

######################################################## MAIN CODE
### prepare annotations - repeatMasker, exons of genes from ENSEMBL annotation, miRBase
## exons
exons_df <-
  exons_gr %>%
  as.data.frame(.) %>%
  tibble::as_tibble(.) %>%
  dplyr::left_join(.,genes_info %>% dplyr::select(gene_id, gene_biotype, gene_name) ,by = "gene_id") %>%
  dplyr::select(seqnames:strand, gene_id, gene_name, gene_biotype) %>%
  tidyr::unite(gene_id, c("gene_name", "gene_id"), sep = "|") %>%
  dplyr::filter(!str_detect(gene_biotype, "miRNA"))

## mature miRNA
mirna_gr <- mirna_gr[mcols(mirna_gr)$type == "miRNA"]
mcols(mirna_gr) <- mcols(mirna_gr)[, c("Name")]
names(mcols(mirna_gr)) <- "gene_id"
mcols(mirna_gr)$gene_biotype <- "miRNA.mature"

## repeatMasker
# add info
rmsk_df %<>% 
  dplyr::mutate(repFamily = stringr::str_replace_na(repFamily),
                gene_id = str_c(repFamily, "|", repName) %>% str_remove("NA\\|"),
                gene_biotype = repClass)

# repeatMasker transposable elements
rmsk_te <-
  rmsk_df %>%
  dplyr::filter(str_detect(repClass, "LINE|SINE|LTR")) %>%
  dplyr::select(-c(repName, repClass, repFamily, rmsk_id)) %>%
  GenomicRanges::GRanges(.)

# repeatMasker other RNA
rmsk_RNA <-
  rmsk_df %>%
  dplyr::filter(str_detect(repClass, "RNA")) %>%
  dplyr::select(-c(repName, repClass, repFamily, rmsk_id)) %>%
  GenomicRanges::GRanges(.)

## ENSEMBL
# ensembl mRNAs
ensembl_mrna <-
  exons_df %>%
  dplyr::filter(gene_biotype == "protein_coding") %>%
  GenomicRanges::GRanges(.)

# # ensembl pseudogene
# ensembl_pseudogene <-
#   exons_df %>%
#   dplyr::filter(str_detect(gene_biotype, "pseudogene")) %>%
#   GenomicRanges::GRanges(.)

# ensembl other RNA
ensembl_RNA <-
  exons_df %>%
  dplyr::filter(!str_detect(gene_biotype, "protein_coding")) %>%
  GenomicRanges::GRanges(.)

# joined other RNA
other_RNA <- c(rmsk_RNA, ensembl_RNA)

# create whole chromosomes as "other" annotation
bam_chr <- chrFinder(sample_bam_path)
other_gr <- GenomicRanges::GRanges(seqnames = bam_chr$chr,
                                   ranges = IRanges(start = 1, end = bam_chr$chrlen),
                                   gene_id = rep("other", nrow(bam_chr)),
                                   gene_biotype = rep("other", nrow(bam_chr)))

### add annotation to list
l.annot <- list(miRNA = mirna_gr,
                TE = rmsk_te,
                mRNA = ensembl_mrna,
                # pseudogene = ensembl_pseudogene,
                RNA_other = other_RNA,
                other = other_gr)

# save RDS
saveRDS(l.annot, file = file.path(outpath, "l.annot.RDS"))
