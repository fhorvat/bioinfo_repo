### INFO: get expression in ENCODE mouse dataset
### DATE: Sun Jun 24 16:14:35 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/ENCODE_2014_Nature_GSE49417/Analysis")

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

library(plotly)
library(htmlwidgets)
library(pcaExplorer)
library(pheatmap)
library(RColorBrewer)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "mutate_cond.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# summarizedOverlaps path
se_path <- list.files(inpath, "*.se.RDS", full.names = T)

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = "ensembl.89.*UCSCseqnames.reducedExons.RDS$", full.names = T)

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = "ensembl.89.*UCSCseqnames.geneInfo.csv$", full.names = T)

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/ENCODE_2014_Nature_GSE49417/Data/Mapped/STAR_mm10"

# stats path
stats_path <- list.files(mapped_path, "log.*stats_and_tracks.csv", full.names = T)

######################################################## READ DATA
# read summarizedExperiment from RDS file
se <- readRDS(file = se_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read info about chosen genes
genes_info <- readr::read_csv(file = genes_info_path)

# read stats
stats_df <- readr::read_csv(file = stats_path)

######################################################## MAIN CODE
# construct sample table
sample_table <- 
  tibble(bam_path = list.files(path = mapped_path, pattern = "*.total.bam$", full.names = T)) %>%
  dplyr::mutate(sample_id = str_remove(basename(bam_path), ".total.bam")) %>%
  dplyr::mutate(tissue = str_remove(sample_id, "s_|") %>% str_remove(., "_.*"),
                age = str_remove_all(sample_id, str_c(tissue, "|s_|r[0-9]{1,}.PE$|_"))) %>%
  dplyr::left_join(stats_df %>% dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA), by = "sample_id") %>%
  dplyr::select(sample_id, tissue, age, library_size, bam_path) %>% 
  dplyr::mutate(tissue = ifelse(tissue == "cns", str_c(tissue, ".", age), tissue)) %T>%
  readr::write_csv(., path = file.path(outpath, "ENCODE_2014_Nature_GSE49417.sample_table.csv"))

# get total length of all exons for each transcript
exons_width <-
  width(exons_gr) %>%
  sum(.) %>%
  tibble(gene_id = names(.), width = .)

# set tissue order
tissue_order <- c("cns.E11.5", "cns.E14", "cns.E18", "frontallobe", "cortex", "cerebellum",
                  "stomach", "liver", "duodenum", "smintestine", "lgintestine", "colon", 
                  "lung", "heart", "bladder", "kidney", "thymus", "mammarygland", "spleen", 
                  "ovary", "testis", "placenta")

### FPKM
# get data.frame of counts, transform to FPKM, write
fpkm_df_avg <-
  assay(se) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  as.tibble(.) %>%
  set_colnames(., str_replace(colnames(.), ".total.bam", "")) %>%
  tidyr::gather(key = sample_id, value = counts, -gene_id) %>%
  dplyr::left_join(., sample_table %>% dplyr::select(sample_id, library_size, tissue), by = "sample_id") %>%
  dplyr::left_join(., exons_width, by = "gene_id") %>%
  dplyr::mutate(library_size = round(library_size / 1E6, 6),
                width = round(width / 1E3, 3),
                fpm = (counts / library_size),
                fpkm = (fpm / width)) %>%
  dplyr::select(gene_id, sample_id, tissue, fpkm) %>%
  dplyr::group_by(gene_id, tissue) %>%
  dplyr::summarise(average_fpkm = mean(fpkm)) %>%
  dplyr::ungroup(.) %>%
  tidyr::spread(key = tissue, value = average_fpkm) %>% 
  dplyr::select_(.dots = (c("gene_id", tissue_order))) %T>%
  readr::write_csv(., path = file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".FPKM.csv")))

### PCA and full heatmap
# get gene_id of protein coding genes
protein_genes <- 
  genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# prepare DDS version of sample table
sample_table_dds <- 
  sample_table %>% 
  dplyr::select(-bam_path, -library_size) %>% 
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)

# filter summarizedExperiment to include only chosen stage and protein coding genes, set colData
se_filt <- se[rownames(se) %in% protein_genes, ]
colData(se_filt) <- DataFrame(sample_table_dds)
se_filt$tissue <- factor(se_filt$tissue)

# make DESeqDataSet
dds <- DESeqDataSet(se_filt, design = ~tissue)

# data for PCA = rlog transformed counts
rlog_df <-
  vst(dds, blind = T) %>%
  assay(.)

# calculate pca
pca <-
  rlog_df %>%
  t(.) %>%
  stats::prcomp(.)

# gets percent of variance for each principal component
percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)

# makes data.frame for ggplot
pca_plot_df <-
  tibble(PC1 = pca$x[, 1],
         PC2 = pca$x[, 2],
         sample_id = colnames(rlog_df)) %>%
  dplyr::left_join(sample_table_dds , by = "sample_id") %>%
  dplyr::mutate(sample_id = str_replace(sample_id, "_r", " ") %>% 
                  str_remove_all(., "s_|.PE|_adult8wks")) %>% 
  dplyr::mutate(tissue = factor(tissue, levels = tissue_order))

# plot PCA
pca_plot <- 
  ggplot(data = pca_plot_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = tissue), size = 7.5, pch = 21, color = "black") +
  # ggrepel::geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  # guides(fill = FALSE, shape = FALSE) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = - 0.2),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  ggsave(filename = file.path(outpath, "PCA.rlog.ENCODE.GRCm38.89.png"), width = 15, height = 15)


### 3D PCA plot
# get colors from 2D plot
colors_pallete <- 
  ggplot_build(pca_plot) %$%
  data[[1]] %>% 
  dplyr::select(fill, group) %>% 
  unique(.) %>% 
  arrange(group) %$%
  fill
  
# data for plot
pca_3Dplot <-
  tibble(PC1 = pca$x[, 1],
         PC2 = pca$x[, 2],
         PC3 = pca$x[, 3],
         sample_id = colnames(rlog_df))  %>%
  dplyr::left_join(sample_table_dds , by = "sample_id") %>%
  dplyr::mutate(sample_id = str_replace(sample_id, "_r", " ") %>% 
                  str_remove_all(., "s_|.PE|_adult8wks")) %>% 
  dplyr::mutate(tissue = factor(tissue, levels = tissue_order))

# make interactive 3D PCA plot
p <-
  plotly::plot_ly(data = pca_3Dplot,
                  x = ~PC1,
                  y = ~PC2,
                  z = ~PC3,
                  color = ~tissue, 
                  colors = colors_pallete) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = str_c("PC1: ", round(percentVar[1] * 100), "% variance")),
                      yaxis = list(title = str_c("PC2: ", round(percentVar[2] * 100), "% variance")), 
                      zaxis = list(title = str_c("PC3: ", round(percentVar[3] * 100), "% variance"))))

# save as html widget
htmlwidgets::saveWidget(plotly::as_widget(p),
                        file = file.path(outpath, "PCA3D.rlog.ENCODE.GRCm38.89.html"),
                        selfcontained = T)


### get loadings - top 10 genes which contribute to PC1 the most
# create annotation data.frame
annotation_df <- 
  genes_info %>% 
  dplyr::select(gene_id, gene_name) %>% 
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "gene_id")

# plot top 10 loadings
par(mar = c(1, 1, 1, 1))
png(file = file.path(outpath, "PCA.rlog.top10_loadings.ENCODE.GRCm38.89.png"), width = 800, height = 600, units = "px")
hi_loadings(pca, topN = 10, annotation = annotation_df)
dev.off()


### distance heatmap
# calculate distance
dist_df <-
  rlog_df %>%
  t(.) %>%
  dist(.)

# make matrix
dist_matrix <- as.matrix(dist_df)
colnames(dist_matrix) <- NULL
rownames(dist_matrix) <- 
  rownames(dist_matrix) %>% 
  str_replace(., "_r", " ") %>% 
  str_remove_all(., "s_|.PE|_adult8wks")

# plot
pheatmap::pheatmap(dist_matrix,
                   clustering_distance_rows = dist_df,
                   clustering_distance_cols = dist_df,
                   col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
                   filename = file.path(outpath, "distHeatmap.rlog.ENCODE.GRCm38.89.png"),
                   height = 15,
                   width = 17)

# ### gene expression heatmap
# # plot
# pheatmap::pheatmap(mat = rlog_df,
#                    cluster_rows = T,
#                    cluster_cols = F,
#                    show_rownames = F, 
#                    show_colnames = T, 
#                    fontsize_col = 20,
#                    col = colorRampPalette(brewer.pal(9, "Blues"))(20),
#                    # breaks = seq(0, 100, by = 10),
#                    filename = file.path(outpath, "expHeatmap.vst.ENCODE.GRCm38.89.png"),
#                    height = 30,
#                    width = 10)
