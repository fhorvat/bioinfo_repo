#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: get FPKM and count tables
### DATE: Thu Apr 25 16:16:56 2019
### AUTHOR: Filip Horvat

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
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# list of read stats and mapped logs
se_path <- args$se_path
sample_table_path <- args$sample_table_path
gtf_path <- args$gtf_path
genes_info_path <- args$genes_info_path
grouping_variables <- args$grouping_variables

######################################################## READ DATA
# summarizedExperiment 
se <- readRDS(se_path)

# read sample table
sample_table <- data.table::fread(sample_table_path)

# read .gtf
gtf_tb <- read_delim(file = gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# read genes info
genes_info <- data.table::fread(genes_info_path)

######################################################## MAIN CODE
### get reduced exons coordinates, calculate length of exons
# convert GTF to GRanges, get only exons, reduce
exons_gr <-
  gtfToGRanges(gtf_tb, filter = "exon") %>%
  GenomicRanges::split(., .$gene_id) %>%
  GenomicRanges::reduce(., ignore.strand = F)

# get total length of all exons for each gene
exons_table <-
  width(exons_gr) %>%
  sum(.) %>%
  tibble(gene_id = names(.), width = .) %>% 
  as.data.table(.)


### calculate FPKMs
# long FPKM/FPM/count table
fpkm_long_tb <-
  assay(se) %>%
  as.data.table(., keep.rownames = "gene_id") %>% 
  melt(., 
       id.vars = c("gene_id"),
       variable.name = "sample_id", 
       value.name = "counts") %>% 
  .[, sample_id := str_remove(sample_id, "\\.bam$")] %>% 
  .[sample_table, on = "sample_id"] %>% 
  .[exons_table, on = "gene_id"] %>%
  .[, fpm := (counts / round(library_size / 1E6, 6))] %>% 
  .[, fpkm := (fpm / round(width / 1E3, 3))] %>% 
  .[, c("gene_id", "sample_id", "counts", "fpm", "fpkm", grouping_variables), with = F]

# save long FPKM table
fpkm_long_tb %>% 
  .[, -grouping_variables, with = F] %>% 
  readr::write_csv(., file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".FPKM_long.csv")))

# wide and save
fpkm_wide_tb <-     
  fpkm_long_tb %>%
  dcast(., gene_id ~ sample_id, value.var = "fpkm") %>% 
  .[genes_info, on = "gene_id", `:=`(gene_name = i.gene_name, 
                                     coordinates_mm10 = str_c(i.seqnames, i.start, i.end, sep = " "), 
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
                                     coordinates_mm10 = str_c(i.seqnames, i.start, i.end, sep = " "), 
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
                                     coordinates_mm10 = str_c(i.seqnames, i.start, i.end, sep = " "), 
                                     strand = i.strand, 
                                     gene_biotype = i.gene_biotype, 
                                     gene_description = i.gene_description)] %>% 
  .[] %T>%
  write_csv(., file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".FPKM_stats.csv")))

