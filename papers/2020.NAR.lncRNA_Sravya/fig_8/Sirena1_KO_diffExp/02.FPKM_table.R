#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: 
### DATE: Thu Apr 25 16:16:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# calculate standard error of the mean FPKM value
standard.error <- function(x) {
  sqrt(var(x) / length(x))
}

######################################################## PATH VARIABLES
### experiment
# set experiment name
experiment <- "Lnc1_KO"
experiment_name <- "Lnc1_KO"

# set grouping variable(s)
grouping_variables <- c("stage", "genotype")


### working dir
# set working directory 
setwd(file.path("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/diffExp/lnc1_KO"))


### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### other experiment paths
# set base experiment path
base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec")

# mapped path
mapped_path <- file.path(base_path, "Data/Mapped/STAR_mm10")

# documentation path
documentation_path <- file.path(base_path, "Data/Documentation")

# analysis path
analysis_path <- inpath


### documentation
# set ensembl version
ensembl_version <- 93

# sample table path
sample_table_path <- file.path(documentation_path, "lnc1_KO.RNAseq.20181211.sampleTable.clean.csv")

# stats and tracks path
stats_and_tracks_path <- list.files(mapped_path, ".*\\.stats_and_tracks\\.csv", full.names = T)

# summarizedExperiment path
se_path <- list.files(path = analysis_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.*.se.RDS$"), full.names = T)


### genome
# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)

######################################################## READ DATA
# read sample table
sample_table <- data.table::fread(sample_table_path)

# read stats and tracks table
stats_and_tracks <- data.table::fread(stats_and_tracks_path)

# summarizedExperiment 
se <- readRDS(se_path)

# read genes info
genes_info <- data.table::fread(genes_info_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

######################################################## MAIN CODE
# join sample table with stats and tracks
sample_table[stats_and_tracks, on = "sample_id", `:=`(library_size = genome.mapped_minus_rDNA)]

# get total length of all exons for each gene
exons_table <-
  width(exons_gr) %>%
  sum(.) %>%
  tibble(gene_id = names(.), width = .) %>% 
  as.data.table(.)

### save count table for GEO submission
# get assay from SummarizedExperiment, tidy column names and save
counts_tb <- 
  assay(se) %>%
  as.data.table(., keep.rownames = "ENSEMBL.93.gene_id") %>% 
  data.table::setnames(., old = 1:ncol(.), 
                       new = colnames(.) %>% 
                         str_remove_all(., "^s_|\\.SE\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam")) %>% 
  readr::write_csv(., file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".raw_counts.csv")))

### calculate FPKMs
# long FPKM/FPM/count table
fpkm_long_tb <-
  assay(se) %>%
  as.data.table(., keep.rownames = "gene_id") %>% 
  melt(., 
       id.vars = c("gene_id"),
       variable.name = "sample_id", 
       value.name = "counts") %>% 
  .[, sample_id := str_remove(sample_id, "\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam|\\.total\\.bam|\\.bam")] %>% 
  .[sample_table, on = "sample_id"] %>% 
  .[exons_table, on = "gene_id"] %>%
  .[, fpm := (counts / round(library_size / 1E6, 6))] %>% 
  .[, fpkm := (fpm / round(width / 1E3, 3))] %>% 
  .[, c("gene_id", "sample_id", "counts", "fpm", "fpkm", grouping_variables), with = F]

# save long FPKM table
fpkm_long_tb %>% 
  .[, -grouping_variables, with = F] %>% 
  readr::write_csv(., file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".FPKM_long.csv")))

# save wide FPKM table
fpkm_wide_tb <-     
  fpkm_long_tb %>%
  dcast(., gene_id ~ sample_id, value.var = "fpkm") %>% 
  .[genes_info, on = "gene_id", `:=`(gene_name = i.gene_name, 
                                     coordinates = str_c(i.seqnames, i.start, i.end, sep = " "), 
                                     strand = i.strand, 
                                     gene_biotype = i.gene_biotype, 
                                     gene_description = i.gene_description)] %>% 
  .[] %T>%
  write_csv(., file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".FPKM.csv")))

# average FPKM across stages and/or genotype 
fpkm_mean_tb <-
  fpkm_long_tb %>%
  .[, .(fpkm_mean = round(mean(fpkm), 3)), by = c("gene_id", grouping_variables)] %>%
  .[, grouped_variables := do.call(str_c, c(.SD, sep = "_")), .SDcols = grouping_variables] %>% 
  dcast(., gene_id ~ grouped_variables, value.var = "fpkm_mean") %>% 
  .[genes_info, on = "gene_id", `:=`(gene_name = i.gene_name, 
                                     coordinates = str_c(i.seqnames, i.start, i.end, sep = " "), 
                                     strand = i.strand, 
                                     gene_biotype = i.gene_biotype, 
                                     gene_description = i.gene_description)] %>% 
  .[] %T>%
  write_csv(., file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".FPKM_mean.csv")))

# long FPKM with mean, SD and standard error for each gene
fpkm_stats <-
  fpkm_long_tb %>%
  .[, list(avg_fpkm = mean(fpkm),
           SD = sd(fpkm),
           SE = standard.error(fpkm)),
    by = c("gene_id", grouping_variables)] %>% 
  .[genes_info, on = "gene_id", `:=`(gene_name = i.gene_name, 
                                     coordinates = str_c(i.seqnames, i.start, i.end, sep = " "), 
                                     strand = i.strand, 
                                     gene_biotype = i.gene_biotype, 
                                     gene_description = i.gene_description)] %>% 
  .[] %T>%
  write_csv(., file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".FPKM_stats.csv")))


### get ratio between KO/WT for chosen genes
# SD(X/Y) = sqrt(SD(X)^2 / Y^2 + ( SD(Y)^2 * X^2 ) / Y^4)
# get flanking genes ID's
flank_genes <- c("ENSMUSG00000030854", "ENSMUSG00000010307", "ENSMUSG00000049516", 
                 "ENSMUSG00000110001", "ENSMUSG00000043262", "ENSMUSG00000014402", 
                 "ENSMUSG00000030851", "ENSMUSG00000063229")

# filter long table - GV
fpkm_flank_genes_GV <- 
  fpkm_long_tb[gene_id %in% flank_genes, 
               list(mean_fpkm = mean(fpkm),SD = sd(fpkm)),
               by = c("gene_id", grouping_variables)] %>% 
  .[stage == "GV", ] %>%
  .[, stage := NULL] %>% 
  dcast(., gene_id ~ genotype, value.var = c("mean_fpkm", "SD")) %>% 
  .[, `:=`(Null_WT_ratio_of_means = (mean_fpkm_Null / mean_fpkm_WT), 
           Null_WT_ratio_of_means_SD = sqrt((SD_Null^2 / mean_fpkm_WT^2) + ((SD_WT^2 * mean_fpkm_Null^2) /  mean_fpkm_WT^4)))] %>% 
  .[] %T>% 
  readr::write_csv(., "flanking_genes.GV.FPKM_ratio.csv")
  
# filter long table - MII
fpkm_flank_genes_MII <- 
  fpkm_long_tb[gene_id %in% flank_genes, 
               list(mean_fpkm = mean(fpkm), SD = sd(fpkm)),
               by = c("gene_id", grouping_variables)] %>% 
  .[stage == "MII", ] %>%
  .[, stage := NULL] %>% 
  dcast(., gene_id ~ genotype, value.var = c("mean_fpkm", "SD")) %>% 
  .[, `:=`(Null_WT_ratio_of_means = (mean_fpkm_Null / mean_fpkm_WT), 
           Null_WT_ratio_of_means_SD = sqrt((SD_Null^2 / mean_fpkm_WT^2) + ((SD_WT^2 * mean_fpkm_Null^2) /  mean_fpkm_WT^4)))] %>% 
  .[] %T>% 
  readr::write_csv(., "flanking_genes.MII.FPKM_ratio.csv")

