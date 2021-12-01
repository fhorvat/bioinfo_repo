### INFO: 
### DATE: Mon Oct 28 18:43:00 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/PhD/algorithms_and_programming/2019_10_28/homework")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# LTR coordinates path
ltr_path <- file.path(inpath, "rmsk.mm10.20180919.MT2A_B.fa.out.csv")

# exon coordinates path
exons_path <- file.path(inpath, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.chr1_all_exons.csv")

######################################################## READ DATA
# read LTRs 
ltr_tb <- readr::read_csv("rmsk.mm10.20180919.MT2A_B.fa.out.csv")

# read exons
exons_tb <- readr::read_csv("ensembl.93.GRCm38.p6.20180919.UCSCseqnames.chr1_all_exons.csv")

######################################################## MAIN CODE
### a)
# transform LTRs table to GRanges
ltr_gr <- 
  ltr_tb %>% 
  GenomicRanges::makeGRangesFromDataFrame(df = ., keep.extra.columns = T, 
                                          seqnames.field = "chr", start.field = "ltr_start", 
                                          end.field = "ltr_end", strand.field = "ltr_strand")
ltr_gr <- ltr_gr[seqnames(ltr_gr) == "chr1"]
seqlevels(ltr_gr) <- "chr1"
mcols(ltr_gr)$unique_ID <- str_c("LTR_", 1:length(ltr_gr))

# transform exon table to GRanges
exons_gr <- 
  exons_tb %>% 
  GenomicRanges::GRanges(.) 

# transform exon regions to GRangesList, split to GRangesList and find ranges of every gene, unlist
gene_gr <- 
  exons_gr %>% 
  split(., mcols(.)$gene_id) %>% 
  range(.) %>% 
  unlist(.)

# find all genes which have at least one LTR insertion in gene body
genes_with_ltrs <- subsetByOverlaps(gene_gr, ltr_gr, ignore.strand = T)

# find all LTRs inserted into genes
ltrs_in_genes <- subsetByOverlaps(ltr_gr, gene_gr, ignore.strand = T)
  

### b)
# find all overlaps between LTRs and genes
overlaps_all <- findOverlaps(ltr_gr, gene_gr, ignore.strand = T)

# get strands of LTRs and genes
ltr_strands <- ltr_gr[queryHits(overlaps_all)] %>% strand(.) %>% as.character(.)
gene_strands <- gene_gr[subjectHits(overlaps_all)] %>% strand(.) %>% as.character(.)

# get number of sense/antisense insertions
sum(ltr_strands == gene_strands)
sum(ltr_strands != gene_strands)


### c)
# get reduced exonic regions as GRanges
exons_reduced <- GenomicRanges::reduce(exons_gr, ignore.strand = T)

# for each gene, get reduced exonic regions
exons_grlist <- 
  split(exons_gr, exons_gr$gene_id) %>% 
  GenomicRanges::reduce(., ignore.strand = F)


### d)
# check if exons GRangesList and genes GRanges objects are parallel
all(names(exons_by_gene_reduced) == names(gene_gr)) 

# find difference between gene coordinates and exon coordinates to get intronic coordinates. Remove empty elements (=monoexonic genes)
introns_grlist <- psetdiff(gene_gr, exons_grlist, ignore.strand = F)
introns_grlist <- introns_grlist[elementNROWS(introns_grlist) > 0]

# subset by overlaps all introns with insertion
introns_with_ltrs <- subsetByOverlaps(introns_grlist, ltr_gr, ignore.strand = T)
length(introns_with_ltrs)


### e)
# One thing to have in mind here that the start of the GRanges object always has to be lower than the end. 
# In other words, ranges can't be negative. This means that start/end of the range cooresponds to 5'/3' end of the biological 
# feature which they represents only if the feature is on the "+" strand. If the feature is on the "-" strand, end of the 
# range cooresponds to 5' end of the biological feature and start of the range to the 3' end. 

# get all LTRs not inserted in the genes
ltrs_intergenic <- ltr_gr[!(mcols(ltr_gr)$unique_ID %in% mcols(ltrs_in_genes)$unique_ID)]

# expand their region 120kb downstream - in the case of + strand this is end + 120 kb. In the case of - strand this is start - 120kb
ltrs_intergenic_downstream_plus <- 
  ltrs_intergenic[strand(ltrs_intergenic) == "+"] %>% 
  resize(., width = 1, fix = "end") %>% 
  GenomicRanges::shift(., shift = 1) %>% 
  resize(., width = 120000, fix = "start")

ltrs_intergenic_downstream_minus <- 
  ltrs_intergenic[strand(ltrs_intergenic) == "-"] %>% 
  resize(., width = 1, fix = "end") %>% 
  GenomicRanges::shift(., shift = -1) %>% 
  resize(., width = 120000, fix = "start")

# join together, subset by overlaps
ltrs_intergenic_downstream <- c(ltrs_intergenic_downstream_plus, ltrs_intergenic_downstream_minus)
subsetByOverlaps(ltrs_intergenic_downstream, gene_gr, ignore.strand = T) %>% length(.)

### f) 
# get the distances to nearest gene on any strand
distanceToNearest(ltrs_intergenic, gene_gr, ignore.strand = T) %>% 
  mcols(.) %$% 
  distance %>% 
  mean(.)

