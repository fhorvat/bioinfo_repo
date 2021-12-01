### INFO: expression of lnc1 in ENCODE data set 
### DATE: Tue Aug 21 14:25:50 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/lncRNA_expression/lnc1_expression")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

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
exons_path <- list.files(path = genome_path, pattern = "ensembl.89.*UCSCseqnames.reducedExons.RDS$", full.names = T)

# gene info path
gene_info_path <- file.path(genome_path, "ensembl.89.GRCm38.p5.20180615.UCSCseqnames.geneInfo.csv")

# path to ENCODE mapped samples
samples_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/ENCODE_2014_Nature_GSE49417/Analysis"

# sample table path
sample_table_path <- list.files(samples_path, pattern = ".*.sample_table.csv", full.names = T)

# summarizedOverlaps path
se_path <- file.path(samples_path, "ensembl.89.GRCm38.p5.20180615.UCSCseqnames.ENCODE_2014_mouse.se.RDS")

######################################################## READ DATA
# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read gene info
gene_info <- readr::read_csv(gene_info_path)

# read sample table
sample_table <- readr::read_csv(sample_table_path)

# reads summarizedOverlaps .RDS
se <- readRDS(file = se_path)

######################################################## MAIN CODE
# lnc1 ENSEMBL ID
lnc1_ensembl <- "ENSMUSG00000110001"

# set tissue order
tissue_order <- c("cns.E11.5", "cns.E14", "cns.E18", "frontallobe", "cortex", "cerebellum",
                  "stomach", "liver", "duodenum", "smintestine", "lgintestine", "colon", 
                  "lung", "heart", "bladder", "kidney", "thymus", "mammarygland", "spleen", 
                  "ovary", "testis", "placenta")

# tissue order with organ systems
tissue_order_df <- 
  tibble(tissue = tissue_order, 
         organ_system = c(rep("nervous", 6), 
                          rep("digestive", 6), 
                          rep("other", 7), 
                          rep("sex", 3)))

# get total length of all exons for each transcript
exons_width <-
  width(exons_gr) %>%
  sum(.) %>%
  tibble(gene_id = names(.), width = .)

# get data.frame of counts, transform to FPKM
fpkm_df <-
  assay(se) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  dplyr::filter(gene_id == lnc1_ensembl) %>% 
  as.tibble(.) %>%
  set_colnames(., str_replace(colnames(.), ".total.bam", "")) %>%
  tidyr::gather(key = sample_id, value = counts, -gene_id) %>%
  dplyr::left_join(., sample_table %>% dplyr::select(sample_id, library_size, tissue), by = "sample_id") %>%
  dplyr::left_join(., exons_width, by = "gene_id") %>%
  dplyr::mutate(library_size = round(library_size / 1E6, 6),
                width = round(width / 1E3, 3),
                fpm = (counts / library_size),
                fpkm = (fpm / width)) %>%
  # dplyr::filter(tissue != "ovary") %>%
  dplyr::select(tissue, fpkm) %>% 
  dplyr::mutate(tissue = factor(tissue, levels = tissue_order))
  
### visualize
# calculate statistics for error bars
stats_df <- 
  fpkm_df %>% 
  dplyr::group_by(tissue) %>% 
  dplyr::summarise(N = n(), 
                   mean = mean(fpkm), 
                   median = median(fpkm),
                   sd = sd(fpkm), 
                   stand_err = sd / sqrt(N)) %>% 
  dplyr::mutate(confidence_interval = stand_err * qt((0.95 / 2) + 0.5, N - 1)) %>% 
  dplyr::left_join(., tissue_order_df, by = "tissue") %>% 
  dplyr::mutate(tissue = factor(tissue, levels = tissue_order))

## BARS
# mean + sd
ggplot(stats_df, aes(x = tissue, y = mean, fill = tissue)) + 
  geom_errorbar(aes(ymin = 0, ymax = mean + sd), width = 0.2) +
  geom_bar(stat = "identity") +
  # coord_cartesian(ylim = c(0, 0.1)) +
  ylab("mean FPKM") + 
  guides(fill = FALSE) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 12, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  ggsave(filename = file.path(outpath, "encode.barplot.mean_sd.nolim.png"), width = 12, height = 10)

# mean + sd flipped
ggplot(stats_df, aes(x = tissue, y = mean, fill = tissue)) + 
  geom_errorbar(aes(ymin = 0, ymax = mean + sd), width = 0.2) +
  geom_bar(stat = "identity") +
  coord_flip(ylim = c(0, 0.1)) +
  ylab("mean FPKM") + 
  guides(fill = FALSE) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 12, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  ggsave(filename = file.path(outpath, "encode.barplot.flip.mean_sd.png"), width = 12, height = 10)

# ## DOTS
# # dotplot - mean + sd
# ggplot() +
#   geom_jitter(data = fpkm_df, aes(x = tissue, y = fpkm, colour = tissue), size = 3) + 
#   stat_summary(data = fpkm_df, aes(x = tissue, y = fpkm),
#                fun.y = mean, fun.ymin = mean, fun.ymax = mean,
#                geom = "crossbar", width = 0.3) +
#   geom_errorbar(data = stats_df, aes(x = tissue, ymin = mean - sd, ymax = mean + sd), width = 0.2) +
#   coord_cartesian(ylim = c(0, 0.1)) +
#   ylab("FPKM") + 
#   guides(color = FALSE) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   ggsave(filename = file.path(outpath, "dots_data.mean_sd.png"), width = 15, height = 10)

# # mean + stand_err
# ggplot(stats_df, aes(x = tissue, y = mean, colour = tissue)) +
#   geom_errorbar(aes(ymin = mean - stand_err, ymax = mean + stand_err), width = 0.1) +
#   geom_point(size = 5, shape = 21, fill = "white") +
#   coord_cartesian(ylim = c(0, 0.1)) +
#   ylab("mean FPKM") +
#   guides(color = FALSE) +
#   theme_bw() +
#   theme(axis.title.x = element_text(size = 10),
#         axis.title.y = element_text(size = 10),
#         axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5),
#         axis.text.y = element_text(size = 12, vjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   ggsave(filename = file.path(outpath, "dots.mean_se.png"), width = 12, height = 10)

# # mean + se
# ggplot(stats_df, aes(x = tissue, y = mean, fill = tissue)) + 
#   geom_bar(stat = "identity") +
#   geom_errorbar(aes(ymin = mean - stand_err, ymax = mean + stand_err), width = 0.2) +
#   # coord_cartesian(ylim = c(0, 0.1)) +
#   ylab("mean FPKM") + 
#   guides(fill = FALSE) +
#   theme_bw() +
#   theme(axis.title.x = element_text(size = 10), 
#         axis.title.y = element_text(size = 10), 
#         axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5), 
#         axis.text.y = element_text(size = 12, vjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) + 
#   ggsave(filename = file.path(outpath, "barplot.mean_se.ovaries.png"), width = 12, height = 10)
