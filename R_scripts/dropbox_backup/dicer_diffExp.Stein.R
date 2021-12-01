#!/home/students/fhorvat/bin/Rscript
### RUN: qsub -q MASTER -M fihorvat@gmail.com -m n -N pbs.Rscript.summarizeOverlaps -l select=ncpus=10:mem=100g -j oe 03f.Dang2016.FPKM.R
### INFO: 
### DATE: Tue May 15 12:58:41 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/Dicer_expression")

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
library(biomaRt)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "gtfToGRanges.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set output path
outpath <- getwd()

# reference path
ref_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# reduced exons path
gtf_path <- file.path(ref_path, "ensembl.91.GRCm38.p5.20180512.UCSCseqnames.gtf.gz")

# gene info path
info_path <- file.path(ref_path, "ensembl.91.GRCm38.p5.20180512.geneInfo.csv")

# stats and tracks, bam paths
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Stein_2015_PLoSGenet_GSE57514/Data/Mapped/STAR_mm10"

######################################################## READ DATA
# read gtf
gtf_df <- read_delim(file = gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# read info about genes
genes_info <- readr::read_csv(info_path)

# get sample table
sample_table <- 
  tibble(bam_path = list.files(mapped_path, pattern = "*.Aligned.sortedByCoord.out.bam$", full.names = T)) %>% 
  dplyr::mutate(ID = basename(bam_path) %>% stringr::str_remove(., ".Aligned.sortedByCoord.out.bam"), 
                genotype = str_remove_all(ID, "s_|_r[123]{1}")) %>% 
  dplyr::filter(str_detect(genotype, "Dicer"))

######################################################## MAIN CODE
# ### get info about genes from Biomart
# # load Mart of mouse database from ensembl
# # sometimes function useMart isn't able to connect to server and returns error, this chunck repeats useMart until there is no error
# mart <- "error"
# count <- 0
# while(class(mart) == "character"){
#   count <- count + 1
#   print(str_c("mmusculus_gene_ensembl", " ", count))
#   mart <- tryCatch(expr = useMart(host = "dec2017.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl"),
#                    error = function(x) return("error"))
# }
# 
# # get info about genes
# ensembl_gene_info <- getBM(attributes = c("ensembl_exon_id", "is_constitutive", "ensembl_gene_id"), mart = mart)

# prepare sample table for dds
sample_table_dds <- 
  sample_table %>% 
  dplyr::select(-bam_path) %>% 
  as.data.frame(.) %>%
  set_rownames(., .$ID)

# get exon ranges, take only protein coding genes
exons_gr <- 
  gtfToGRanges(gtf_df, filter = "exon")  %>% 
  asdf(.) %>% 
  as.tibble(.) %>% 
  dplyr::filter(gene_biotype == "protein_coding")

# get Dicer1 gene ID
dicer_id <-
  genes_info %>%
  dplyr::filter(gene_name == "Dicer1") %$%
  gene_id

# get IDs of first 10 exons in Dicer1
dicer_10exons <- 
  exons_gr %>% 
  dplyr::filter(transcript_id == "ENSMUST00000041987") %>% 
  dplyr::arrange(desc(start)) %>% 
  dplyr::slice(1:10) %$% 
  exon_id

# remove all exons except first 10 in Dicer1, reduce ranges
exons_gr_filt <- 
  exons_gr %>% 
  dplyr::filter(!((gene_id == dicer_id) & !(exon_id %in% dicer_10exons))) %>% 
  GenomicRanges::GRanges(.) %>% 
  GenomicRanges::split(., .$gene_id) %>%
  GenomicRanges::reduce(., ignore.strand = T)

# ### counts
# # get count of reads, save summarizedExperiment as RDS
# bamfiles <- Rsamtools::BamFileList(sample_table$bam_path, yieldSize = 2000000)
# BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))
# se <- GenomicAlignments::summarizeOverlaps(features = exons_gr_filt,
#                                            reads = bamfiles,
#                                            mode = "Union",
#                                            singleEnd = TRUE,
#                                            ignore.strand = TRUE)
# saveRDS(se, file = file.path(outpath, "ensembl.GRCm38.91.Dicer1.exons1_to_10.summarizedOverlaps.RDS"))

# read summarizedExperiment from RDS file
se <- readRDS(file = file.path(outpath, "ensembl.GRCm38.91.Dicer1.exons1_to_10.summarizedOverlaps.RDS"))
colnames(se) <- str_remove(colnames(se), ".Aligned.sortedByCoord.out.bam")
colData(se) <- DataFrame(sample_table_dds)
se$genotype <- factor(se$genotype, levels = c("Dicer_WT", "Dicer_KO"))

# make DESeqDataSet
dds <- DESeqDataSet(se, design = ~genotype)

# run DESeq
dds_deseq <- DESeq(dds)

# get results, shrink logFC
dds_shrink <- lfcShrink(dds_deseq, contrast = c("genotype", "Dicer_KO", "Dicer_WT"))

# get Dicer results
dicer_res <- 
  dds_shrink %>% 
  as.data.frame(.) %>% 
  tibble::rownames_to_column(., var = "gene_id") %>% 
  as.tibble(.) %>% 
  dplyr::filter(gene_id == dicer_id) %T>%
  readr::write_csv(x = ., path = file.path(outpath, "diffExp.Dicer1.exons1_to_10.Stein_DicerKOvsWT.csv"))

