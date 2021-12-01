#!/common/WORK/fhorvat/programi/R/R-3.4.3/bin/Rscript
### RUN: qsub -q MASTER -M fihorvat@gmail.com -m n -N pbs.Rscript.summarizeOverlaps -l select=ncpus=10:mem=100g -j oe 03e.Xue2013.FPKM.R
### INFO: get expression of all genes in Xue data
### DATE: Mon Mar 12 14:51:03 2018
### AUTHOR: Filip Horvat

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

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set output path
outpath <- getwd()

# reduced exons path
exons_path <- "/common/WORK/fhorvat/reference/human/hg38/ensembl.GRCh38.91.20180312.UCSCnames.clean.reducedExons.RDS"

# stats and tracks, bam paths
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Xue_2013_Nature_GSE44183/Data/Mapped/STAR_hg38"

# info about genes path
genes_info_path <- "/common/WORK/fhorvat/reference/human/hg38/ensembl.GRCh38.91.20180312.UCSCnames.clean.geneInfo.csv"

# info about CNOT path 
CNOT_info_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis/Mus_musculus.GRCm38.89.20180305.CNOTinfo.csv"

######################################################## READ DATA
# read ENSEMBL reduced exons 
exons_gr <- readRDS(file = exons_path)

# read stats and tracks table
stats_df <- readr::read_csv(file = list.files(mapped_path, pattern = "*stats_and_tracks.csv", full.names = T))

# read additional info about genes
genes_info <- readr::read_csv(genes_info_path)

# read additional info about CNOT
CNOT_info <- readr::read_csv(CNOT_info_path)

######################################################## MAIN CODE
# get total length of all exons for each transcript
exons_width <- 
  width(exons_gr) %>%
  sum(.) %>% 
  tibble(gene_id = names(.), width = .)

# construct sample table
sample_table <- 
  tibble(bam_path = list.files(path = mapped_path, pattern = "*.bam$", full.names = T)) %>% 
  dplyr::mutate(sample_id = str_replace(basename(bam_path), ".total.bam", "")) %>% 
  dplyr::left_join(stats_df %>% dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA), by = "sample_id")

# get gene_id of protein coding genes
protein_genes <- 
  genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# # get count of reads, save summarizedExperiment as RDS
# bamfiles <- Rsamtools::BamFileList(sample_table$bam_path, yieldSize = 2000000)
# BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))
# se <- GenomicAlignments::summarizeOverlaps(features = exons_gr,
#                                            reads = bamfiles,
#                                            mode = "Union",
#                                            singleEnd = FALSE,
#                                            ignore.strand = TRUE)
# saveRDS(se, file = file.path(outpath, basename(exons_path) %>% str_replace(., ".reducedExons.RDS", ".Xue2013.se.RDS")))

# read summarizedExperiment from RDS file
se <- readRDS(file = file.path(outpath, "ensembl.GRCh38.91.20180312.UCSCnames.clean.Xue2013.se.RDS")) 


# ### FPKM - average across stage
# # get data.frame of counts, transform to FPKM
# fpkm_avg <-
#   assay(se) %>%
#   as.data.frame(.) %>%
#   tibble::rownames_to_column(., var = "gene_id") %>%
#   as.tibble(.) %>%
#   set_colnames(., str_replace(colnames(.), ".total.bam", "")) %>%
#   tidyr::gather(key = sample_id, value = counts, -gene_id) %>%
#   dplyr::left_join(., sample_table %>% dplyr::select(sample_id, library_size), by = "sample_id") %>%
#   dplyr::left_join(., exons_width, by = "gene_id") %>%
#   dplyr::mutate(library_size = round(library_size / 1E6, 6),
#                 width = round(width / 1E3, 3),
#                 fpm = (counts / library_size),
#                 fpkm = (fpm / width)) %>%
#   dplyr::select(gene_id, sample_id, fpkm) %>%
#   dplyr::mutate(sample_id = str_replace(sample_id, ".PE", ""),
#                 stage = str_replace_all(sample_id, "^s_|_r[0-9]{1,}", "")) %>%
#   dplyr::group_by(gene_id, stage) %>%
#   dplyr::summarise(average_fpkm = mean(fpkm)) %>%
#   dplyr::ungroup(.) %>%
#   tidyr::spread(key = stage, value = average_fpkm) %T>%
#   readr::write_csv(., path = file.path(outpath, "ensembl.GRCh38.91.20180312.UCSCnames.clean.Xue2013.avgFPKM.csv"))
# 
# # get average FPKM of CNOT complex genes
# fpkm_avg_filt <-
#   fpkm_avg %>%
#   dplyr::right_join(., CNOT_info %>% dplyr::select(gene_id = hsapiens_homolog_ensembl_gene, gene_name), by = "gene_id") %>%
#   dplyr::select(gene_name, everything(), -gene_id) %>%
#   dplyr::arrange(gene_name) %T>%
#   readr::write_csv(., path = file.path(outpath, "CNOT.GRCh38.91.20180312.UCSCnames.clean.Xue2013.avgFPKM.csv"))


### FPKM - individual samples
# get data.frame of counts, transform to FPKM
fpkm_ind <-
  assay(se) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  as.tibble(.) %>%
  set_colnames(., str_replace(colnames(.), ".total.bam", "")) %>%
  tidyr::gather(key = sample_id, value = counts, -gene_id) %>%
  dplyr::left_join(., sample_table %>% dplyr::select(sample_id, library_size), by = "sample_id") %>%
  dplyr::left_join(., exons_width, by = "gene_id") %>%
  dplyr::mutate(library_size = round(library_size / 1E6, 6),
                width = round(width / 1E3, 3),
                fpm = (counts / library_size),
                fpkm = (fpm / width)) %>%
  dplyr::select(gene_id, sample_id, fpkm) %>%
  dplyr::mutate(sample_id = str_replace(sample_id, ".PE", "")) %>%
  tidyr::spread(key = sample_id, value = fpkm)

# # get individual FPKM of CNOT complex genes
# fpkm_ind_filt <-
#   fpkm_ind %>%
#   dplyr::right_join(., CNOT_info %>% dplyr::select(gene_id = hsapiens_homolog_ensembl_gene, gene_name), by = "gene_id") %>%
#   dplyr::select(gene_name, everything(), -gene_id) %>%
#   dplyr::arrange(gene_name) %T>%
#   readr::write_csv(., path = file.path(outpath, "CNOT.GRCh38.91.20180312.UCSCnames.clean.Xue2013.FPKM.csv"))
# 
# # plot CNOT FPKM values across samples as boxplot
# fpkm_ind_filt %>%
#   tidyr::gather(key = sample_id, value = counts, -gene_name) %>%
#   dplyr::mutate(stage = str_replace_all(sample_id, "^s_|_r[0-9]{1,}", "")) %>%
#   dplyr::filter(str_detect(gene_name, "Cnot"),
#                 !str_detect(stage, "8C_embryo|pronucleus")) %>%
#   dplyr::mutate(gene_name = factor(gene_name, levels = str_c("Cnot", c(1:6, "6l", 7:11))),
#                 stage = factor(stage, levels = c("MII", "1C", "2C_blastomere", "4C_blastomere", "8C_blastomere", "morula"))) %>%
#   ggplot(data = ., aes(x = stage, counts)) +
#   geom_boxplot() +
#   geom_point(aes(colour = stage),  size = 2) +
#   scale_x_discrete(drop = FALSE, position = "bottom") +
#   facet_grid(gene_name ~ ., scales = "free") +
#   theme_bw() +
#   theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks.x = element_blank()) +
#   ggsave(filename = file.path(outpath, "results", "boxplot.CNOT.GRCh38.91.20180312.UCSCnames.clean.Xue2013.FPKM.png"), height = 20, width = 10, limitsize = F)


### PCA
# # prepare sample table for dds
# sample_table_dds <- 
#   sample_table %>% 
#   dplyr::select(-bam_path, -library_size) %>% 
#   dplyr::mutate(sample_id = str_replace(sample_id, ".PE", ""), 
#                 stage = str_replace_all(sample_id, "^s_|_r[0-9]{1,}", "")) %>%
#   as.data.frame(.) %>%
#   magrittr::set_rownames(., .$sample_id)
#
# ## rlog transformed countS PCA
# # filter summarizedExperiment to include only chosen stage and protein coding genes, set colData
# se_filt <- se[rownames(se) %in% protein_genes, ]
# colData(se_filt) <- DataFrame(sample_table_dds)
# 
# # make DESeqDataSet
# dds <- DESeqDataSet(se_filt, design = ~stage)
# 
# # data for PCA
# PCA_df <-
#   rlog(dds, blind = T) %>%
#   assay(.)
# 
# # calculates pca
# pca <-
#   PCA_df %>%
#   t(.) %>%
#   stats::prcomp(.)
# 
# # gets percent of variance for each principal component
# percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)
# 
# # makes data.frame for ggplot, plots PCA
# tibble(PC1 = pca$x[, 1],
#        PC2 = pca$x[, 2],
#        sample_id = colnames(PCA_df)) %>%
#   dplyr::left_join(sample_table_dds , by = "sample_id") %>%
#   dplyr::mutate(sample_id = str_replace(sample_id, "s_|r", "") %>% str_replace_all(., "_", " ")) %>%
#   ggplot(data = ., aes(x = PC1, y = PC2, label = sample_id, color = stage)) +
#   geom_point(aes(color = stage), size = 5) +
#   # geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
#   xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
#   ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
#   guides(color = guide_legend(override.aes = list(size = 5))) +
#   theme_bw() +
#   theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
#         axis.title.y = element_text(size = 15, vjust = - 0.2),
#         axis.text.x = element_text(size = 15),
#         axis.text.y = element_text(size = 15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# save plot
# ggsave(filename = file.path(outpath, "results", str_c("PCAplot.Xue2013.GRCh38.91.rlog.png")), 
#        plot = PCA_plot, width = 10, height = 10)

# ## individual FPKM PCA
# # data for PCA
# PCA_df <-
#   fpkm_ind %>%
#   dplyr::filter(gene_id %in% protein_genes) %>%
#   dplyr::mutate_all(.funs = funs(log2(. + 1))) %>%
#   as.data.frame(.) %>%
#   tibble::column_to_rownames(var = "gene_id")
# 
# # calculates pca
# pca <-
#   PCA_df %>%
#   t(.) %>%
#   stats::prcomp(.)
# 
# # gets percent of variance for each principal component
# percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)
# 
# # makes data.frame for ggplot, plots PCA
# PCA_plot <-
#   tibble(PC1 = pca$x[, 1],
#          PC2 = pca$x[, 2],
#          sample_id = colnames(PCA_df)) %>%
#   dplyr::mutate(stage = str_replace_all(sample_id, "^s_|_r[0-9]{1,}", ""),
#                 sample_id = str_replace(sample_id, "s_|r", "") %>% str_replace_all(., "_", " ")) %>%
#   ggplot(data = ., aes(x = PC1, y = PC2, label = sample_id, color = stage)) +
#   geom_point(aes(color = stage), size = 5) +
#   # geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
#   xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
#   ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
#   guides(color = guide_legend(override.aes = list(size = 5))) +
#   theme_bw() +
#   theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
#         axis.title.y = element_text(size = 15, vjust = - 0.2),
#         axis.text.x = element_text(size = 15),
#         axis.text.y = element_text(size = 15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# # save plot
# ggsave(filename = file.path(outpath, "results", str_c("PCAplot.Xue2013.GRCh38.91.FPKM.png")), 
#        plot = PCA_plot, width = 10, height = 10)

## separate genes on PCA
# # data for PCA
# PCA_df <-
#   fpkm_ind %>% 
#   dplyr::filter(gene_id %in% protein_genes) %>% 
#   dplyr::mutate_at(.vars = vars(-gene_id), .funs = funs(log2(. + 1))) %>% 
#   as.data.frame(.) %>% 
#   tibble::column_to_rownames(var = "gene_id")


# ### t-sne
# library(Rtsne)
# 
# # prepare data
# tsne_df <- 
#   fpkm_ind %>% 
#   dplyr::filter(gene_id %in% protein_genes) %>%
#   .[, 2:ncol(.)] %>% 
#   unique() %>% 
#   t(.)
# 
# ## Executing the algorithm on curated data
# tsne <- Rtsne(tsne_df, dims = 2, perplexity = 3, verbose = TRUE, max_iter = 500)
# 
# # makes data.frame for ggplot, plots PCA
# tsne_plot <-
#   tibble(x = tsne$Y[, 1],
#          y = tsne$Y[, 2],
#          sample_id = rownames(tsne_df)) %>%
#   dplyr::mutate(stage = str_replace_all(sample_id, "^s_|_r[0-9]{1,}", ""),
#                 sample_id = str_replace(sample_id, "s_|r", "") %>% str_replace_all(., "_", " ")) %>%
#   ggplot(data = ., aes(x = x, y = y, label = sample_id, color = stage)) +
#   geom_point(aes(color = stage), size = 5) +
#   # geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
#   guides(color = guide_legend(override.aes = list(size = 5))) +
#   theme_bw() +
#   theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
#         axis.title.y = element_text(size = 15, vjust = - 0.2),
#         axis.text.x = element_text(size = 15),
#         axis.text.y = element_text(size = 15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# # save plot
# ggsave(filename = file.path(outpath, "results", str_c("tsne.Xue2013.GRCh38.91.FPKM.png")),
#        plot = tsne_plot, width = 10, height = 10)
