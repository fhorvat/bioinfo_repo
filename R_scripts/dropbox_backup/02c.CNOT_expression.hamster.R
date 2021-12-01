#!/common/WORK/fhorvat/programi/R/R-3.4.3/bin/Rscript
### INFO: get expression of all CNOT genes in Fugaku's data
### DATE: Fri Mar 09 22:34:51 2018
### AUTHOR: Filip Horvat
### RUN: qsub -q MASTER -M fihorvat@gmail.com -m n -N pbs.Rscript.summarizeOverlaps -l select=ncpus=10:mem=300g -j oe 01.expression_analysis.R

rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)

######################################################## PATH VARIABLES
# set output path
outpath <- getwd()

# reduced exons path
exons_path <- "/common/WORK/fhorvat/reference/hamster_golden/mesAur1/ensembl/ensembl.MesAur1.0.91.20180309.reducedExons.RDS"

# info about genes path
genes_info_path <- file.path(outpath, "Mus_musculus.GRCm38.89.20180305.CNOTinfo.csv")

# stats and tracks, bam paths
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mesAur1_ensembl"

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## READ DATA
# read ENSEMBL reduced exons 
exons_gr <- readRDS(file = exons_path)

# read additional info about genes
genes_info <- readr::read_csv(genes_info_path)

# read stats and tracks table
stats_df <- readr::read_csv(file = list.files(mapped_path, pattern = "*stats_and_tracks.csv", full.names = T))

######################################################## MAIN CODE
# filter exons by genes_info table
exons_gr <- exons_gr[names(exons_gr) %in% genes_info$btaurus_gene_id]

# get total length of all exons for each transcript
exons_width <- 
  width(exons_gr) %>%
  sum(.) %>% 
  tibble(gene_id = names(.), width = .)

# construct sample table
sample_table <- 
  tibble(bam_path = list.files(path = mapped_path, pattern = "*.bam$", full.names = T)) %>% 
  dplyr::mutate(sample_id = str_replace(basename(bam_path), ".genome.Aligned.sortedByCoord.out.bam", "")) %>% 
  dplyr::left_join(stats_df %>% dplyr::select(sample_id, genome.mapped_minus_rDNA), by = "sample_id") %>% 
  dplyr::rename(library_size = genome.mapped_minus_rDNA)

# get count of reads, save summarizedExperiment as RDS
bamfiles <- Rsamtools::BamFileList(sample_table$bam_path, yieldSize = 2000000)
BiocParallel::register(BiocParallel::MulticoreParam())
se <- GenomicAlignments::summarizeOverlaps(features = exons_gr,
                                           reads = bamfiles,
                                           mode = "Union",
                                           singleEnd = TRUE,
                                           ignore.strand = TRUE)
saveRDS(se, file = file.path(outpath, "Graf_2014.UMD3.1.91.bosTau8.reducedExons.summarizedOverlaps.RDS"))

# read summarizedExperiment from RDS file
se <- readRDS(file = file.path(outpath, "Graf_2014.UMD3.1.91.bosTau8.reducedExons.summarizedOverlaps.RDS")) 

### FPKM
# get data.frame of counts, transform to FPKM
fpkm_df <-
  assay(se) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  as.tibble(.) %>%
  set_colnames(., str_replace(colnames(.), ".genome.Aligned.sortedByCoord.out.bam", "")) %>%
  tidyr::gather(key = sample_id, value = counts, -gene_id) %>%
  dplyr::left_join(., sample_table %>% dplyr::select(sample_id, library_size), by = "sample_id") %>%
  dplyr::left_join(., exons_width, by = "gene_id") %>%
  dplyr::mutate(library_size = round(library_size / 1E6, 6),
                width = round(width / 1E3, 3),
                fpm = (counts / library_size),
                fpkm = (fpm / width)) %>%
  dplyr::select(gene_id, sample_id, fpkm) %>%
  dplyr::mutate(stage = str_replace_all(sample_id, "s_|_GSxB.*|_GS.*", "")) %>% 
  dplyr::group_by(gene_id, stage) %>% 
  dplyr::summarise(average_fpkm = mean(fpkm)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::spread(key = stage, value = average_fpkm) %>%
  dplyr::left_join(., genes_info %>% dplyr::select(gene_id = btaurus_gene_id, gene_name), by = "gene_id") %>% 
  dplyr::select(gene_name, everything()) %>% 
  dplyr::select(-gene_id) %>% 
  dplyr::arrange(gene_name) %T>%
  readr::write_csv(., path = file.path(outpath, "FPKM.Graf_2014.CNOT.ensembl.UMD3.1.91.bosTau8.csv"))

# plot
fpkm_df %>%
  tidyr::gather(key = sample_id, value = fpkm, -gene_name) %>%
  ggplot(data = ., aes(x = gene_name, y = fpkm, fill = gene_name)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_grid(sample_id ~ .) +
  ggsave(file.path(outpath, "test.pdf"), width = 10, height = 30)
