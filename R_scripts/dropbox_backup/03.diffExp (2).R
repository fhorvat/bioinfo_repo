### INFO: lnc1 KO differential expression analysis
### DATE: Sun Mar 11 04:54:44 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/expression")

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
library(ggrepel)
library(RColorBrewer)
library(plotly)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path 
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# reduced exons path
exons_path <- list.files(path = genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.reducedExons.RDS$", full.names = T)

# gene info path
gene_info_path <- list.files(genome_path, "ensembl.91.GRCm38.p5.*.UCSCseqnames.geneInfo.csv", full.names = T)

# sample table path
sample_table_path <- list.files(inpath, pattern = ".*filtered.sample_table.csv", full.names = T)

# mean genotype fpkm path
fpkm_path <- file.path(inpath, "lncRNA_KO.mean_genotype.FPKM.csv")
  
# summarizedOverlaps path
se_path <- file.path(inpath, "ensembl.91.GRCm38.p5.20180512.lncKO.summarizedOveralaps.RDS")

######################################################## READ DATA
# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read info about all ENSEMBL annotated genes
gene_info <- 
  readr::read_csv(gene_info_path) %>% 
  tidyr::unite(cooridnates, seqnames, start,  end, sep = " ") 

# sample table path
sample_table <- readr::read_csv(file = sample_table_path)

# read FPKM table
fpkm_df <- readr::read_csv(fpkm_path)

# read summarizedExperiment from RDS file
se <- readRDS(file = se_path)

######################################################## MAIN CODE
### prepare data
# prepare sample table for DESeq colData
sample_table_dds <- 
  sample_table %>% 
  dplyr::select(-bam_path, -library_size) %>% 
  # dplyr::mutate(genotype = ifelse(genotype == "WT", str_c(genotype, "_", str_extract(sample_id, "201[6,7]{1}[:alpha:]{3}")), genotype)) %>% 
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)

# get gene_id of protein coding genes
protein_genes <- 
  gene_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# get total length of all exons for each transcript
exons_width <-
  width(exons_gr) %>%
  sum(.) %>%
  tibble(gene_id = names(.), width = .)

# filter summarizedExperiment to include only chosen stage and protein coding genes, set colData
se_filt <- se[rownames(se) %in% protein_genes, ]
se_filt <- se_filt[, match(rownames(sample_table_dds), colnames(se_filt))]
colData(se_filt) <- DataFrame(sample_table_dds)

# make DESeqDataSet
dds <- DESeqDataSet(se_filt, design = ~genotype)

# ### PCA plot
# # data for PCA = rlog transformed counts
# rlog_df <-
#   rlog(dds, blind = T) %>%
#   assay(.)
# 
# # calculates pca
# pca <-
#   rlog_df %>%
#   t(.) %>%
#   stats::prcomp(.)
# 
# # gets percent of variance for each principal component
# percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)
# 
# # makes data.frame for ggplot, plots PCA
# pca_plot <-
#   tibble(PC1 = pca$x[, 1],
#          PC2 = pca$x[, 2],
#          sample_id = colnames(rlog_df)) %>%
#   dplyr::left_join(sample_table_dds , by = "sample_id") %>%
#   dplyr::mutate(sample_id = str_remove_all(sample_id, "s_|.reseq|") %>% str_replace_all(., "_", " ") %>% str_replace(., ".SE.", " ")) %>%
#   ggplot(data = ., aes(x = PC1, y = PC2, label = sample_id, color = genotype)) +
#   geom_point(aes(color = genotype), size = 7.5) +
#   # geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
#   xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
#   ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
#   theme_bw() +
#   theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
#         axis.title.y = element_text(size = 15, vjust = - 0.2),
#         axis.text.x = element_text(size = 15),
#         axis.text.y = element_text(size = 15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) + 
#   ggsave(filename = file.path(outpath, "PCAplot.lncKO.PC1_2.rlog.png"), width = 10, height = 10)

### DESeq2
# run DESeq
dds_deseq <- DESeq(dds)

#### lnc1Null vs. WT
## get results
# set results to fetch
comparison <- c("Lnc1Null", "WT")

# get results, shrink logFC
dds_shrink <- lfcShrink(dds_deseq, contrast = c("genotype", comparison))

## MA plot
# get table for plot
plot_df <-
  dds_shrink %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  tibble::as.tibble(.) %>% 
  dplyr::left_join(., gene_info %>% dplyr::select(gene_id, gene_name), by = "gene_id") %>% 
  dplyr::arrange(padj) %>%
  dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_name) %>% 
  dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
                sign = ifelse(padj < 0.1, "yes", "no"),
                regulation = ifelse(lfc > 0, "up", "down"), 
                regulation = replace(regulation, sign == "no", "not_sign"), 
                regulation = factor(regulation, levels = c("not_sign", "up", "down"))) %>% 
  dplyr::arrange(regulation)

# get labels
labels_df <- 
  plot_df %>% 
  dplyr::filter(sign == "yes")

# plot
ggplot() + 
  geom_point(data = plot_df, aes(x = mean, y = lfc, color = regulation), size = 3, shape = 20) +
  geom_label_repel(data = labels_df, aes(x = mean, y = lfc, label = gene_name), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
  scale_x_log10(limits = c(1e-01, 1e5), 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(limits = c(-7, 5),
                     breaks = c(-7:5)) +
  scale_colour_manual(values = c(not_sign = "gray30", up = "red3", down = "blue3")) +
  guides(color = FALSE) +
  ggtitle("lnc1NULL vs WT") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  ggsave(filename = file.path(outpath, "MAplot.lnc1NUll_vs_WTall.png"), width = 10, height = 10)



