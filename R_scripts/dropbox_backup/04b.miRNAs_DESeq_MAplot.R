### INFO: 
### DATE: Tue Jan 07 15:36:14 2020
### AUTHOR: vfranke, modified by Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2019_Dec/Analysis/expression/miRBase")

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
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)
library(DESeq2)
library(ggrepel)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gtf path
gtf_path <- file.path(genome_path, "miRBase.22.mm10.20181605.gff3")

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2019_Dec"

# FPM path
fpm_path <- file.path(inpath, "miRBase.22.mm10.20181605.FPM.csv")

# sample table path
sample_table_path <- file.path(base_path, "Data/Documentation")
sample_table_path <- list.files(sample_table_path, ".*\\.sampleTable\\.csv", full.names = T)

# library size path
library_size_path <- file.path(base_path, "Data/Mapped/STAR_mm10/7_perfect_reads") 
library_size_path <- list.files(library_size_path, "library_sizes.txt", full.names = T)

# pepa's annoatation path
pepa_path <- file.path(inpath, "mirAnnot.dt.rda")

# featureCounts path
counts_path <- file.path(inpath, "miRBase.22.mm10.20181605.counts.txt")

######################################################## READ DATA
# read gtf info
mirbase_gr <- rtracklayer::import(gtf_path)

# read small RNA-seq expression
fpm_tb <- readr::read_csv(fpm_path)

# read sample table
sample_tb <- readr::read_csv(sample_table_path)

# load Pepa's annotation
load(pepa_path) 

# read counts from featureCounts
counts_tb <- readr::read_delim(counts_path, delim = "\t", comment = "#") 

######################################################## MAIN CODE
### prepare tables
# get miRBase annotation table
mirbase_tb <- 
  mirbase_gr %>% 
  as_tibble(.) %>% 
  tidyr::unite(coordinates, seqnames, start, end, sep = " ") %>% 
  dplyr::select(gene_id = Name, ID = ID, parent = Derives_from, type, coordinates) %>% 
  dplyr::filter(type == "miRNA") %>% 
  dplyr::left_join(., mirAnnot.dt %>% dplyr::select(ID, strand_type = type), by = "ID") %>% 
  dplyr::mutate(arm = str_extract(gene_id, "3p|5p"),
                arm = replace(arm, is.na(arm), "not_defined"))

# get counts with ID's
counts_tb %<>%
  set_colnames(., basename(colnames(.))) %>% 
  tidyr::unite(coordinates, Chr, Start, End, sep = " ", remove = F) %>% 
  dplyr::left_join(., mirbase_tb %>% dplyr::select(ID, coordinates), by = "coordinates")

# prepare sample table for DESeq colData
sample_table_dds <-
  sample_tb %>%
  # dplyr::filter(genotype != "HET") %>% 
  # dplyr::filter(!(str_detect(sample_id, "KO_7B_r4|KO_2_r5"))) %>% 
  dplyr::select(sample_id, genotype) %>% 
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)


## summarized experiment
# counts table
se <-
  counts_tb %>%
  dplyr::select(-c(Geneid:Length)) %>%
  dplyr::select(gene_id = ID, everything()) %>%
  dplyr::mutate_if(is.numeric, round, digits = 0) %>%
  as.data.frame(.) %>%
  set_rownames(., .$gene_id) %>%
  dplyr::select(-gene_id) %>%
  as.matrix(.)

# filter summarizedExperiment to include only chosen stage, set colData
se_filt <- se
se_filt <- SummarizedExperiment(list(counts = se_filt))
colnames(se_filt) <- str_remove_all(colnames(se_filt), "\\.24to31nt|\\.21to23nt|\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam|\\.total\\.bam|\\.bam$")
se_filt <- se_filt[, colnames(se_filt)[match(rownames(sample_table_dds), colnames(se_filt))]]

# check if colnames of assay match rownames of sample table DataFrame
if(all(colnames(se_filt) == rownames(sample_table_dds))){
  
  # set sample table as colData
  colData(se_filt) <- DataFrame(sample_table_dds)
  
}else{
  
  # stop script with warrning
  stop("Columns in assay are not matching row of sample table. Please check your data annotation")
  
}


## DDS
# make DESeqDataSet
dds <- DESeqDataSet(se_filt, design = ~genotype)

# DESeq
dds_deseq <- DESeq(dds)

# get results, shrink logFC
# dds_shrink <- lfcShrink(dds_deseq, contrast = c("genotype", "KO", "WT"))
dds_shrink <- results(dds_deseq, contrast = c("genotype", "KO", "WT"))

# get results table, add to sheet
results_df <-
  dds_shrink %>%
  as_tibble(., rownames = "ID") %>%
  dplyr::arrange(padj) 


### plot as MA plot
# get axis limits
results_limits <-
  results_df %>%
  dplyr::summarise(x_limit = baseMean %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.),
                   y_limit = log2FoldChange %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.))

# data for plot
plot_df <-
  results_df %>%
  dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, ID) %>%
  dplyr::left_join(., mirbase_tb %>% dplyr::select(ID, strand_type, arm), by = "ID") %>% 
  dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                padj = replace(padj, padj == 0, .Machine$double.xmin),
                sign = ifelse(padj < 0.1, "yes", "no"),
                regulation = replace(arm, padj >= 0.1, "not_defined"),
                regulation = factor(regulation, levels = c("not_defined", "3p", "5p"))) %>%
  dplyr::arrange(regulation)

# plot
ma_plot <-
  ggplot() +
  geom_point(data = plot_df, aes(x = mean, y = lfc, color = regulation, alpha = regulation), size = 2.5, shape = 20) +
  scale_x_log10(limits = c(0.01, results_limits$x_limit),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(limits = c(-results_limits$y_limit, results_limits$y_limit),
                     breaks = c(-results_limits$y_limit:results_limits$y_limit)) +
  scale_colour_manual(labels = c("not_defined" = "not significant", "5p" = "5'", "3p" = "3'"),
                      values = c("not_defined" = "gray50", "5p" = "#1a75ff", "3p" = "red2")) +
  scale_alpha_manual(values = c("not_defined" = 0.5, "5p" = 1, "3p" = 1)) +
  guides(alpha = F) +
  xlab("mean expression") +
  ylab(str_c("log2 fold change: ", "Dicer X/X", " vs. ", "WT E15.5", "\n")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)) +
  theme(legend.title = element_blank())

# # turns off axis titles and legend
# ma_plot <-
#   ma_plot +
#   theme(legend.position = "none") +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank())

# save plot
ggsave(filename = file.path(outpath,
                            str_c("miRBase",
                                  "plot", "MA", "DESeq2", "genotype",
                                  str_c("DicerX_KO", "_vs_", "WT"),
                                  "png", sep = ".")),
       plot = ma_plot, width = 12, height = 10)


### add gene IDs as labels
# set p-adjusted and log2FC cutoff
padj_cut <- 0.1
lfc_cut <- 1.5

# filter data
plot_df_labels <-
  plot_df %>%
  dplyr::mutate(padj_sign = ifelse(padj <= padj_cut, T, F),
                lfc_sign = ifelse((abs(lfc) >= lfc_cut), T, F)) %>%
  dplyr::filter(padj_sign & lfc_sign)

# annotation table
annotations <- tibble(xpos = Inf,
                      ypos = -Inf,
                      annotateText = str_c("label cutoff: ",
                                           "p-adjusted <= ", padj_cut,
                                           ", log2FC >= ", lfc_cut))

# add labels
ma_plot_labeled <-
  ma_plot +
  geom_text(data = plot_df_labels,
            aes(x = mean, y = lfc, label = ID),
            check_overlap = TRUE, size = 3, hjust = 0, vjust = 1.5,
            colour = "black", fontface = "plain") +
  geom_text(data = annotations, aes(x = xpos, y = ypos, label = annotateText),
            colour = "black", fontface = "italic", size = 2.5,
            hjust = 1.03, vjust = -0.5)

# save plot
ggsave(filename = file.path(outpath,
                            str_c("miRBase",
                                  "plot", "MA", "DESeq2", "genotype",
                                  "labeled", 
                                  str_c("DicerX_KO", "_vs_", "WT"),
                                  "png", sep = ".")),
       plot = ma_plot_labeled, width = 12, height = 10)
