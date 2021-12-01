### INFO: 
### DATE: Sat Jun 29 15:02:03 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Fugaku_intron_CpG/CpG_frequency")

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
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- file.path(getwd(), "../count_intronic_reads")

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# FPKM path
fpkm_path <- file.path(inpath, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.intronic.FPKM.csv")

# 1-cell bigWig coverage path
coverage_1C_path <- file.path(inpath, "../intronic_reads/s_1cell.WE.PE.intronic.bw")

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read intronic FPKM
fpkm_tb <- readr::read_csv(fpkm_path)

# read 1C coverage
coverage_1C <- rtracklayer::import.bw(coverage_1C_path)

######################################################## MAIN CODE
### get coordinates of introns where reads are sensitive to DBA
# filter fpkm table
intronic_coords <- 
  fpkm_tb %>% 
  dplyr::filter(s_GV.WE.PE == 0, 
                s_1cell.WE.PE >= 1, 
                s_1cell.WE_DNAm.PE <= s_1cell.WE.PE * 0.1) %>% 
  dplyr::select(gene_id, coordinates, strand, s_1cell.WE.PE) 

# create GRanges
intronic_whole_gr <- 
  intronic_coords %>% 
  tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.)

# set names
names(intronic_whole_gr) <- mcols(intronic_whole_gr)$gene_id
mcols(intronic_whole_gr)$gene_id <- NULL


### coordinates of 1-cell reads
# reduce coverage of 1-cell reads, intersect with whole DBA sensitive introns
intronic_1C_gr <- 
  coverage_1C %>%
  reduce(.) %>% 
  intersect(., intronic_whole_gr, ignore.strand = T)

# findOverlaps with DBA senstive introns
intronic_gr_hits <- 
  findOverlaps(intronic_1C_gr, intronic_whole_gr, ignore.strand = T) %>% 
  subjectHits(.) %>% 
  intronic_whole_gr[.]

# add names and strand information
names(intronic_1C_gr) <- names(intronic_gr_hits)
strand(intronic_1C_gr) <- strand(intronic_gr_hits)


### remove 1C reads covered regions from whole introns
# set difference
intronic_without_1C <- setdiff(intronic_whole_gr, intronic_1C_gr, ignore.strand = F)

# findOverlaps with DBA senstive introns
intronic_gr_hits <- 
  findOverlaps(intronic_without_1C, intronic_whole_gr, ignore.strand = T) %>% 
  subjectHits(.) %>% 
  intronic_whole_gr[.]

# add names
names(intronic_without_1C) <- names(intronic_gr_hits)


### get only introns which are not completely covered in 1C reads
# remove from 1C coverage
intronic_1C_gr <- intronic_1C_gr[names(intronic_1C_gr) %in% names(intronic_without_1C)]


### CpG frequency 
# get frequency of CpG dinucleotide in whole introns
intronic_without_1C_CpG <- 
  intronic_without_1C %>% 
  getSeq(x = Mmusculus, .) %>% 
  as.character(.) %>% 
  tibble(gene_id = names(.), 
         seq = .) %>% 
  dplyr::mutate(cpg_count = str_count(string = seq, pattern = "CG"), 
                width = nchar(seq)) %>% 
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(cpg_count = sum(cpg_count), 
                   width = sum(width)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(cpg_freq = cpg_count / width) %>% 
  dplyr::select(gene_id, cpg_freq_intron = cpg_freq, cpg_count_intron = cpg_count, width_intron = width)

# get frequency of CpG dinucleotide in parts of introns covered by 1C reads
intronic_1C_CpG <- 
  intronic_1C_gr %>% 
  getSeq(x = Mmusculus, .) %>% 
  as.character(.) %>% 
  tibble(gene_id = names(.), 
         seq = .) %>% 
  dplyr::mutate(cpg_count = str_count(string = seq, pattern = "CG"), 
                width = nchar(seq)) %>% 
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(cpg_count = sum(cpg_count), 
                   width = sum(width)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(cpg_freq = cpg_count / width) %>% 
  dplyr::select(gene_id, cpg_freq_1C = cpg_freq, cpg_count_1C = cpg_count, width_1C = width)


### plot
# join intronic CpG with 1C reads CpG
CpG_tb <- 
  left_join(intronic_CpG, intronic_1C_CpG, by = "gene_id") %>% 
  left_join(., intronic_coords, by = "gene_id")

# shape data for plot
CpG_plot_tb <- 
  CpG_tb %>% 
  # dplyr::filter(width_intron > 75,
  #               width_1C > 75) %>%
  dplyr::select(gene_id, cpg_freq_intron, cpg_freq_1C) %>% 
  tidyr::gather(key = category, value = cpg_freq, -gene_id) %>% 
  dplyr::mutate(category = factor(category, levels = c("cpg_freq_intron", "cpg_freq_1C")))

# plot
ggplot(CpG_plot_tb, aes(x = category, y = cpg_freq, color = category)) +
  geom_boxplot(outlier.colour = NULL) +
  geom_point() +
  geom_line(aes(group = gene_id), color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank()) +
  ggsave(file.path(outpath, "CpG_frequency.introns_vs_1C_reads.png"), width = 10, height = 10)

