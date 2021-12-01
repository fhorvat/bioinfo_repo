### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/Papd7_KO/retrotransposon_expression_signature")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

library(GenomicRanges)
library(karyoploteR)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
expandRange = function(x, upstream=2000, downstream=1000) {
  strand_is_minus = strand(x) == "-"
  on_plus = which(!strand_is_minus)
  on_minus = which(strand_is_minus)
  start(x)[on_plus] = start(x)[on_plus] - upstream
  start(x)[on_minus] = start(x)[on_minus] - downstream
  end(x)[on_plus] = end(x)[on_plus] + downstream
  end(x)[on_minus] = end(x)[on_minus] + upstream
  x
}

######################################################## PATH VARIABLES
# inpath
inpath <- getwd()

# outpath
outpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
gene_info_path <- list.files(genome_path, "ensembl\\.93.*\\.UCSCseqnames\\.geneInfo\\.csv", full.names = T)

# repeat masker path
rmsk_path <- list.files(genome_path, "rmsk.*\\.clean\\.fa\\.out\\.gz", full.names = T)

# differential expression results path
diffExp_path <- file.path(inpath, "DE-analysis-20191211.xlsx")

# retrotransposon expression path
expression_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/datasets/2019_Dec.Papd7/Analysis/retrotransposon_expression"
retro_mean_path <- file.path(expression_path, "rmsk.mm10.20180919.joined_rmsk_id.FPKM_mean.csv")

######################################################## READ DATA
# read gene info
gene_info <- readr::read_csv(gene_info_path)

# read repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read differential expression results
diffExp_tb <- openxlsx::read.xlsx(diffExp_path) %>% as_tibble(.)

# read retrotransposon expression
retro_mean_tb <- readr::read_csv(retro_mean_path)

######################################################## MAIN CODE
### prepare and clean the data
# clean the differential expression table
diffExp_tb_tidy <- 
  diffExp_tb %>% 
  dplyr::mutate(gene_id = str_remove(gene_id, "\\.[0-9]+$")) %>% 
  dplyr::left_join(., gene_info %>% dplyr::select(gene_id, seqnames, start, end, strand), by = "gene_id") %>% 
  dplyr::filter(log2FoldChange >= 0)

# create GRanges
genes_gr <- GRanges(diffExp_tb_tidy)

# clean retrotransposon expression
retro_tb_tidy <- 
  retro_mean_tb %>% 
  dplyr::select(gene_id, FPKM_KO = `testes_P7_2-10_KO`, FPKM_WT = testes_WT, 
                coordinates, strand, gene_biotype, gene_description)

# create GRanges
retro_gr <- 
  retro_tb_tidy %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.)


### check how many upregulated genes have upregulated elements near them
# expand ranges of genes
genes_expand_gr <- expandRange(genes_gr, upstream = 10000, downstream = 10000)

# find overlaps with elements
overlaps <- findOverlaps(genes_expand_gr, retro_gr, ignore.strand = T)

# for each gene find expression of the elements nearby
genes_retro_tb <- 
  tibble(gene_id = mcols(genes_expand_gr[queryHits(overlaps)])$gene_id, 
         rmsk_id = mcols(retro_gr[subjectHits(overlaps)])$gene_id, 
         FPKM_KO = mcols(retro_gr[subjectHits(overlaps)])$FPKM_KO, 
         FPKM_WT = mcols(retro_gr[subjectHits(overlaps)])$FPKM_WT) %>% 
  dplyr::mutate(log2FC.KO_vs_WT = log2(FPKM_KO / FPKM_WT)) %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::filter(FPKM_KO >= 1) %>%
  dplyr::filter(log2FC.KO_vs_WT > 1)

# get genes which fit the filtering criteria
genes_filt <- 
  diffExp_tb_tidy %>% 
  dplyr::filter(gene_id %in% genes_retro_tb$gene_id)


### get the same thing from random genes
# get all genes as GRanges
genes_all_gr <- GRanges(gene_info)

# repeat 1000 times
set.seed(1234)
length_random <- purrr::map(1:1000, function(n){
  
  # sample random genes
  random_gr <- genes_all_gr[sample(1:length(genes_all_gr), nrow(diffExp_tb_tidy), replace = F)]
  
  # expand 
  random_expand_gr <- expandRange(random_gr, upstream = 10000, downstream = 10000)
  
  # find overlaps with elements
  overlaps <- suppressWarnings(findOverlaps(random_expand_gr, retro_gr, ignore.strand = T))
  
  # for each gene find expression of the elements nearby
  genes_retro_tb <- 
    tibble(gene_id = mcols(random_expand_gr[queryHits(overlaps)])$gene_id, 
           rmsk_id = mcols(retro_gr[subjectHits(overlaps)])$gene_id, 
           FPKM_KO = mcols(retro_gr[subjectHits(overlaps)])$FPKM_KO, 
           FPKM_WT = mcols(retro_gr[subjectHits(overlaps)])$FPKM_WT) %>% 
    dplyr::mutate(log2FC.KO_vs_WT = log2(FPKM_KO / FPKM_WT)) %>% 
    dplyr::group_by(gene_id) %>% 
    dplyr::filter(FPKM_KO >= 1) %>%
    dplyr::filter(log2FC.KO_vs_WT > 1)
  
  # get mean distance
  return(length(unique(genes_retro_tb$gene_id)))
  
})

# find quantiles
number_tb <- 
  tibble(number = c(nrow(genes_filt), 
                    unlist(length_random)), 
         source = c("our_result", 
                    rep("random", length(length_random)))) %>% 
  dplyr::mutate(promile = ntile(number, 1000)) %>% 
  dplyr::mutate(percentile = ntile(number, 100)) %>% 
  dplyr::arrange(number) %>% 
  dplyr::filter(source == "our_result")


### visualize on chromosomes
# prepare data
subset <- list(1:7, 8:14, 15:21)
chrs <- c(str_c("chr", 1:19), "chrX", "chrY")

# create table for karotype
gene_karyo <-
  diffExp_tb_tidy %>%
  dplyr::select(chr = seqnames, start, end, strand) %>%
  dplyr::mutate(y = runif(nrow(.), -0.75, 0.75)) %>%
  dplyr::group_by(chr) %>%
  split(., .$strand)

### plot dots
for(n in 1:length(subset)){
  
  ### open file
  png(file.path(outpath, str_c("karyo.testis.Papd7_2_10_P14_KO_vs_WT", 
                               "diffExp_upregulated_genes", n, ".png")), 
      width = 1000, height = 1000, units = "px")
  
  # create basic plot
  kp <- plotKaryotype(genome = "mm10", plot.type = 2, chromosomes = chrs[subset[[n]]], cex = 3)
  
  # add plus strand
  kpDataBackground(kp, r0 = 0, r1 = 0.75, color = "#e6e6e6", data.panel = 1)
  kpPoints(kp, chr = gene_karyo[[1]]$chr, x = gene_karyo[[1]]$start, y = gene_karyo[[1]]$y,
           ymin = -1, ymax = 1, col = "red", pch = ".", cex = 5,
           r0 = 0, r1 = 0.75, data.panel = 1)
  
  # add minus strand
  kpDataBackground(kp, r0 = 0, r1 = 0.75, color = "#e6e6e6", data.panel = 2)
  kpPoints(kp, chr = gene_karyo[[2]]$chr, x = gene_karyo[[2]]$start, y = gene_karyo[[2]]$y,
           ymin = -1, ymax = 1, col = "red", pch = ".", cex = 5,
           r0 = 0, r1 = 0.75, data.panel = 2)
  
  ### close file
  dev.off()
}




# ### plot density
# for(n in 1:length(subset)){
#   
#   png(file.path(outpath, str_c("mm10.LINE1_genome_insertions.density.", n, ".png")), width = 1000, height = 1000, units = "px",)
#   kp <- plotKaryotype(genome = "mm10", plot.type = 2, chromosomes = chrs[subset[[n]]], cex = 3)
#   
#   kpPlotDensity(kp, rmsk_tidy[strand(rmsk_tidy) == "+"], 
#                 window.size = 10e5,
#                 r0 = 0, r1 = 0.75, data.panel = 1)
#   kpPlotDensity(kp, exons_gr[strand(exons_gr) == "+"], 
#                 window.size = 10e5,
#                 col="#ddaacc", 
#                 r0 = 0, r1 = 0.75, data.panel = 1)
#   
#   kpPlotDensity(kp, rmsk_tidy[strand(rmsk_tidy) == "-"], 
#                 window.size = 10e5,
#                 r0 = 0, r1 = 0.75, data.panel = 2)
#   kpPlotDensity(kp, exons_gr[strand(exons_gr) == "-"], 
#                 window.size = 10e5,
#                 col="#ddaacc", 
#                 r0 = 0, r1 = 0.75, data.panel = 2)
#   dev.off()
#   
# }
# 
