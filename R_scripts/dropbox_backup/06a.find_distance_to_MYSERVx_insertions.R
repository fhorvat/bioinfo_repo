### INFO: 
### DATE: Fri Dec 11 17:05:38 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/Mov10l1_KO_analysis/testis.RNA_seq/distance_to_retrotransposons")

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
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genomic files
genome_path <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# repeatMasker path
rmsk_path <- file.path(genome_path, "rmsk.Siomi.20200701.clean.fa.out.gz")

# annotation path
annotation_path <- file.path(genome_path, "annotation/Liftoff/MesAur1/ENSEMBL", "ensembl.99.MesAur1.0.20200415.UCSCseqnames.Siomi.liftoff.filt.geneInfo.csv")

######################################################## READ DATA
# read repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read annotation
annotation_tb <- readr::read_csv(annotation_path)

######################################################## MAIN CODE
### tidy data
# filter repeatMasker, get GRanges
rmsk_gr <- 
  rmsk_tb %>% 
  # dplyr::filter(str_detect(repName, "MYSERV|MMERVK")) %>% 
  dplyr::filter(str_detect(repName, "MYSERV")) %>% 
  GRanges(.)

# get gene annotation GRanges
annotation_gr <- 
  annotation_tb %>%
  dplyr::filter(gene_biotype == "protein_coding") %>% 
  GRanges(.)


### read data
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets"

# experiment
experiment_name <- "hamster_testis_Mov10l.8.5dpp.run_2.RNAseq"
experiment_name <- "hamster_testis_Mov10l.0.5dpp.RNAseq"

# results path
results_path <- file.path(base_path, experiment_name,
                          "Analysis/expression.added_PIWIL3.stranded", "results.ensembl.99.MesAur1.0.20200415.UCSCseqnames.Piwil3_Tex101.from_RefSeq")
results_path <- file.path(results_path, "protein_coding.diffExp.DESeq2.genotype_age.significant_results.xlsx")

# read differently expressed genes
results_tb <- openxlsx::read.xlsx(results_path) %>% as_tibble(.)

# get Siomi coordinates of upregulated genes
results_gr <- 
  results_tb %>% 
  dplyr::filter(log2FoldChange > 0) %>% 
  dplyr::select(gene_id, log2FoldChange, baseMean, coordinates) %>%
  dplyr::left_join(., annotation_tb, by = "gene_id") %>% 
  dplyr::filter(!is.na(seqnames)) %>% 
  GRanges(.)

# for each gene find nearest MYSERV/MMERVK insertion
nearest_rmsk <- GenomicRanges::distanceToNearest(results_gr, rmsk_gr, ignore.strand = T)

# create results table
results_distance_tb <- 
  results_gr[queryHits(nearest_rmsk)] %>% 
  as_tibble(.) %>% 
  dplyr::mutate(nearest_insertion_distance = mcols(nearest_rmsk)$distance) %>% 
  tidyr::unite(coordinates.Siomi, seqnames, start, end, sep = " ") %>% 
  dplyr::select(-width) %>% 
  dplyr::mutate(gene_description = str_remove(gene_description, " \\[.*")) %>% 
  dplyr::arrange(nearest_insertion_distance)

# calculate the mean distance 
mean_distance_data <- mean(results_distance_tb$nearest_insertion_distance)
names(mean_distance_data) <- experiment_name

### get distance from random genes
# repeat 1000 times
set.seed(1234)
mean_distances_random <- purrr::map(1:1000, function(n){
  
  # sample random genes
  random_gr <- annotation_gr[sample(1:length(annotation_gr), length(results_gr), replace = F)]
  
  # distance
  nearest_random <- GenomicRanges::distanceToNearest(random_gr, rmsk_gr, ignore.strand = T)
  
  # get mean distance
  return(mean(mcols(nearest_random)$distance))
  
})

# find quantiles
distances_tb <- 
  tibble(distances = c(mean_distance_data, 
                       unlist(mean_distances_random)), 
         source = c(names(mean_distance_data), 
                    rep("random", length(mean_distances_random)))) %>% 
  dplyr::mutate(promile = ntile(distances, 1000)) %>% 
  dplyr::mutate(percentile = ntile(distances, 100)) %>% 
  dplyr::arrange(distances) %>% 
  dplyr::filter(source == names(mean_distance_data))


# save
readr::write_csv(results_insertions_tb, path = file.path(outpath, str_c("oocyte_Mov10l1.KO_vs_WT.upregulated_genes.distance_to_nearest_L1_IAP.csv")))


