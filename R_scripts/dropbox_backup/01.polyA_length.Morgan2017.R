### INFO: comparison of genes polyA tail length and regulation in CNOT6L KO
### DATE: Fri Jun 15 10:27:48 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/review/polyA_tail_length")

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

# ENSEMBL annotated genes info path
ensembl_genes_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.89.GRCm38.p5.20180615.UCSCseqnames.geneInfo.csv"

# CNOT6L FPKM path
CNOT6L_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis/ensembl.GRCm38.89.CNOT6L.avgFPKM.csv"

# results path
results_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis/results"

# polyA lengths paths
polya_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/review/polyA_tail_length/TAIL-seq_Morganetal2017.csv"

######################################################## READ DATA
# read info about all ENSEMBL annotated genes
ensembl_genes_info <- readr::read_csv(ensembl_genes_path)

# read average CNOT6L FPKM
CNOT6L_fpkm <- readr::read_csv(file = CNOT6L_path)

# read CNOT6L significantly diff. exp. genes
CNOT6L_results <- 
  list.files(results_path, "diffExp.CNOT6L.*.signif.csv", full.names = T) %>% 
  purrr::map(., readr::read_csv) %>% 
  dplyr::bind_rows(.)

# read polyA tails lengths
polya_lengths <- readr::read_csv(polya_path)

######################################################## MAIN CODE
# get gene_id of protein coding genes
protein_coding <- 
  ensembl_genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# filter results table
results <- 
  CNOT6L_results %>% 
  dplyr::mutate(stage = str_extract(comparison, "1C|MII|GV"), 
                regulation = ifelse(log2FoldChange > 0, "up", "down")) %>% 
  dplyr::select(gene_id, stage, regulation, logFC_CNOT6L = log2FoldChange, gene_description, GV_KO_FPKM = GV_KO, MII_KO_FPKM = MII_KO)

# filter polyA lengths table, calculate average per gene
polya_df <- 
  polya_lengths %>% 
  dplyr::filter(cell_tissue == "GV", treatment == "CTL") %>% 
  dplyr::select(accession, gene_name = gene, sample, tail_length = poly_a_tail_length) %>%
  dplyr::group_by(accession, gene_name) %>% 
  dplyr::summarise(tail_length = mean(tail_length) %>% round(., 2)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(ensembl_genes_info %>% dplyr::select(gene_name, gene_id), by = "gene_name") %>% 
  dplyr::filter(gene_id %in% protein_coding)

# filter FPKM table, join with polyA lengths
fpkm_df <- 
  CNOT6L_fpkm %>% 
  dplyr::filter(gene_id %in% protein_coding) %>% 
  data.table::setnames(., 2:ncol(.), str_c("s.", colnames(.)[2:ncol(.)])) %>% 
  dplyr::left_join(., polya_df %>% dplyr::select(gene_id, tail_length), by = "gene_id") %>% 
  dplyr::select(gene_id, s.MII_WT, s.MII_KO, tail_length) %>% 
  dplyr::filter(!is.na(tail_length))

# create ratio data.frame
ratio_df <-
  fpkm_df %>% 
  dplyr::select(x = tail_length, y1 = s.MII_KO, y2 = s.MII_WT, gene_id) %>% 
  # dplyr::filter(y2 > 1) %>%
  dplyr::mutate(y = log2(y1 / y2)) %>%
  dplyr::select(gene_id, x, y) %>%
  dplyr::filter(complete.cases(.),
                !is.infinite(x),
                !is.infinite(y)) %>% 
  dplyr::left_join(., results %>% dplyr::filter(stage == "MII") %>% dplyr::select(gene_id, logFC = logFC_CNOT6L, regulation), by = "gene_id") %>%
  dplyr::mutate(regulation = replace(regulation, is.na(regulation), "no"), 
                regulation = factor(regulation, levels = c("no", "down", "up"))) %>%
  dplyr::arrange(regulation) %>% 
  dplyr::filter(regulation == "up")

# plot - blue: #1a75ff
plot_ratio_1 <-
  ratio_df %>%
  ggplot(data = ., aes(x = x, y = y, color = regulation)) +
  geom_point(size = 2.5) +
  scale_x_continuous(limits = c(0, 90), breaks = seq(0, 80, 10)) +
  scale_y_continuous(limits = c(- 8.1, 8.1), breaks = seq(-8, 8, 2)) +
  scale_color_manual(values = c(up = "red3", down = "#1a75ff", no = "gray60"), breaks = c("up", "no")) +
  guides(color = FALSE) +
  geom_hline(yintercept = 0) +
  xlab("polyA tail length (Morgan 2017)") +
  ylab("log2 (MII_KO / MII_WT) CNOT6L") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = 0.3),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c("scatter_ratio.polyA.CNOT6L.MII.KO_vs_WT.noFilter.halfRation.png")), 
       plot = plot_ratio_1, width = 7.5, height = 15)

### histogram 
# all genes data.frame
ratio_df <-
  fpkm_df %>% 
  dplyr::select(tail_length, gene_id) %>% 
  dplyr::left_join(., results %>% dplyr::filter(stage == "MII") %>% dplyr::select(gene_id, logFC = logFC_CNOT6L, regulation), by = "gene_id") %>%
  dplyr::mutate(regulation = replace(regulation, is.na(regulation), "no"), 
                regulation = factor(regulation, levels = c("no", "down", "up"))) %>%
  dplyr::arrange(regulation)

# calculate bin sizes
binwidth <- 10
n <- nrow(ratio_df)
x <- ratio_df[["tail_length"]]
origin <- 0
bin <- (x - origin) %/% binwidth

# add bins
bins_all <- 
  ratio_df %>% 
  dplyr::mutate(bin = bin + 1) %>% 
  dplyr::group_by(bin) %>% 
  dplyr::summarise(size_all = n())

### upregulated genes
# filter upregulated
ratio_df_upregulated <- 
  ratio_df %>% 
  dplyr::filter(regulation == "up")

# calculate bin sizes
binwidth <- 10
n <- nrow(ratio_df)
x <- ratio_df_upregulated[["tail_length"]]
origin <- 0
bin <- (x - origin) %/% binwidth

# add bins
bins_up <- 
  ratio_df_upregulated %>% 
  dplyr::mutate(bin = bin + 1) %>% 
  dplyr::group_by(bin) %>% 
  dplyr::summarise(size_up = n())

# join, calculate fraction
bins_df <- 
  bins_all %>% 
  dplyr::left_join(., bins_up, by = "bin") %>% 
  dplyr::mutate(size_up = replace(size_up, is.na(size_up), 0), 
                fraction = size_up / size_all, 
                bin = bin * 10) %T>% 
  readr::write_csv(., file.path(outpath, str_c("histogram.polyA.CNOT6L.MII.KO_vs_WT.fract.csv")))

# plot
plot_hist <-
  ggplot(data = bins_df, aes(x = bin, y = fraction)) +
  geom_bar(stat='identity') +
  scale_x_continuous(limits = c(0, 90), breaks = seq(0, 80, 10)) +
  guides(color = FALSE) +
  geom_hline(yintercept = 0) +
  xlab("polyA tail length (Morgan 2017)") +
  # ylab("log2 (MII_KO / MII_WT) CNOT6L") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = 0.3),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c("histogram.polyA.CNOT6L.MII.KO_vs_WT.noFilter.png")), 
       plot = plot_hist, width = 15, height = 15)


# ### other options
# # plot 2 - transparent points
# plot_ratio_2 <-
#   ratio_df %>%
#   ggplot(data = ., aes(x = x, y = y)) +
#   geom_point(size = 2.5, alpha = 0.1) +
#   scale_y_continuous(limits = c(- 8.1, 8.1), breaks = seq(-8, 8, 2)) +
#   guides(color = FALSE) +
#   geom_hline(yintercept = 0) +
#   xlab("polyA tail length (Morgan 2017)") +
#   ylab("log2 (MII_KO / GV_KO) CNOT6L") +
#   theme_bw() +
#   theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
#         axis.title.y = element_text(size = 15, vjust = 0.3),
#         axis.text.x = element_text(size = 15),
#         axis.text.y = element_text(size = 15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# # plot 3 - density contours
# plot_ratio_3 <-
#   ratio_df %>%
#   ggplot(data = ., aes(x = x, y = y)) +
#   geom_point(size = 2.5, alpha = 0.1) + 
#   geom_density_2d() + 
#   scale_y_continuous(limits = c(- 8.1, 8.1), breaks = seq(-8, 8, 2)) +
#   guides(color = FALSE) +
#   geom_hline(yintercept = 0) +
#   xlab("polyA tail length (Morgan 2017)") +
#   ylab("log2 (MII_KO / GV_KO) CNOT6L") +
#   theme_bw() +
#   theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
#         axis.title.y = element_text(size = 15, vjust = 0.3),
#         axis.text.x = element_text(size = 15),
#         axis.text.y = element_text(size = 15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# # plot 4 - filled density contours
# plot_ratio_4 <-
#   ratio_df %>%
#   ggplot(data = ., aes(x = x, y = y)) +
#   stat_density_2d(aes(fill = ..level..), geom = 'polygon') +
#   scale_fill_viridis(name = "density") +
#   geom_point(size = 2.5, shape = ".") + 
#   scale_y_continuous(limits = c(- 8.1, 8.1), breaks = seq(-8, 8, 2)) +
#   guides(color = FALSE) +
#   geom_hline(yintercept = 0) +
#   xlab("polyA tail length (Morgan 2017)") +
#   ylab("log2 (MII_KO / GV_KO) CNOT6L") +
#   theme_bw() +
#   theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
#         axis.title.y = element_text(size = 15, vjust = 0.3),
#         axis.text.x = element_text(size = 15),
#         axis.text.y = element_text(size = 15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# # plot 5 - hexbins
# plot_ratio_5 <-
#   ratio_df %>%
#   ggplot(data = ., aes(x = x, y = y)) +
#   geom_hex(bins = 50) +
#   scale_fill_viridis() +
#   geom_point(shape = ".", col = 'white') + 
#   scale_y_continuous(limits = c(- 8.1, 8.1), breaks = seq(-8, 8, 2)) +
#   guides(color = FALSE) +
#   geom_hline(yintercept = 0) +
#   xlab("polyA tail length (Morgan 2017)") +
#   ylab("log2 (MII_KO / GV_KO) CNOT6L") +
#   theme_bw() +
#   theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
#         axis.title.y = element_text(size = 15, vjust = 0.3),
#         axis.text.x = element_text(size = 15),
#         axis.text.y = element_text(size = 15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# # plot 6 - rugs
# plot_ratio_6 <-
#   ratio_df %>%
#   ggplot(data = ., aes(x = x, y = y)) +
#   geom_point(size = 2.5, color = "gray60") +
#   geom_rug(alpha = 0.01) +
#   scale_y_continuous(limits = c(- 8.1, 8.1), breaks = seq(-8, 8, 2)) +
#   guides(color = FALSE) +
#   geom_hline(yintercept = 0) +
#   xlab("polyA tail length (Morgan 2017)") +
#   ylab("log2 (MII_KO / GV_KO) CNOT6L") +
#   theme_bw() +
#   theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
#         axis.title.y = element_text(size = 15, vjust = 0.3),
#         axis.text.x = element_text(size = 15),
#         axis.text.y = element_text(size = 15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# # save plot
# all_plots <- cowplot::plot_grid(plot_ratio_1, plot_ratio_2, plot_ratio_3, plot_ratio_4, plot_ratio_5, plot_ratio_6,
#                                 ncol = 3, labels = "AUTO", align = "v", axis = "lr")
# cowplot::save_plot(file.path(outpath, str_c("test.png")), all_plots, base_height = 15, base_aspect_ratio = 1.5)
