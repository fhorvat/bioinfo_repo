### INFO: plot cross-hair plot of ratio between polyA and total RNA
### DATE: Tue Jul 09 12:41:09 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/polyA_total_ratio")

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

library(plotly)
library(htmlw)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# experiment names
experiment_total_name <- "Fugaku"
experiment_polyA_name <- "Freimer_2018_CurrBiol_GSE92761"

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/"

# get path to FPKM tables
FPKM_total_path <- file.path(base_path, experiment_total_name, "Analysis/expression", 
                             str_c("ensembl.93.GRCm38.p6.20180919.UCSCseqnames.", str_remove(experiment_total_name, "_.*"), ".FPKM_mean.csv"))
FPKM_polyA_path <- file.path(base_path, experiment_polyA_name, "Analysis/expression", 
                             str_c("ensembl.93.GRCm38.p6.20180919.UCSCseqnames.", str_remove(experiment_polyA_name, "_.*"), ".FPKM_mean.csv"))

# path to histones genes
histones_genes_path <- file.path(inpath, "histone_selection.csv")
  
######################################################## READ DATA
# read FPKM tables
FPKM_total <- readr::read_csv(FPKM_total_path)
FPKM_polyA <- readr::read_csv(FPKM_polyA_path)

# read table of histones genes
histones_genes_tb <- readr::read_csv(histones_genes_path)

######################################################## MAIN CODE
# set lnc1 gene
lnc1_gene <- "ENSMUSG00000110001"

# set dormant genes = Mos, Plat, Cyclin B1, Orc6l, and Dcp1a
dormant_genes <- c("ENSMUSG00000078365", "ENSMUSG00000031538", "ENSMUSG00000041431", "ENSMUSG00000031697", "ENSMUSG00000021962")

# get histone genes
histones_gene <- 
  FPKM_total %>%
  dplyr::filter(gene_name %in% histones_genes_tb$gene_name) %$% 
  gene_id


### tidy and join tables
# total table
FPKM_total_tidy <-
  FPKM_total %>% 
  dplyr::select(gene_id, GV.total = GV.WE, MII.total = MII.WE, gene_biotype, gene_name, gene_description)

# polyA table
FPKM_polyA_tidy <-
  FPKM_polyA %>% 
  dplyr::select(gene_id, GV.polyA = GV, MII.polyA = MII, gene_biotype)


### get one table and plot 
# join tables
FPKM_tb <- 
  left_join(FPKM_total_tidy, FPKM_polyA_tidy, by = c("gene_id", "gene_biotype")) %>% 
  dplyr::mutate(gene_biotype = replace(gene_biotype, gene_id == lnc1_gene, "lnc1"), 
                gene_biotype = replace(gene_biotype, gene_id %in% dormant_genes, "dormant"), 
                gene_biotype = replace(gene_biotype, gene_id %in% histones_gene, "histone")) %>% 
  dplyr::filter(gene_biotype %in% c("lnc1", "dormant", "histone", "protein_coding")) %>% 
  dplyr::select(gene_id, GV.polyA, GV.total, MII.polyA, MII.total, gene_biotype, gene_name, gene_description) %>% 
  dplyr::filter_at(.vars = vars(GV.total, MII.total), .vars_predicate = all_vars(. > 1))
  

### construct and save cross-hair plot
# table for plot
plot_tb <- 
  FPKM_tb %>% 
  dplyr::mutate(GV_ratio = log2((GV.polyA + 0.1) / (GV.total + 0.1)), 
                MII_ratio = log2((MII.polyA + 0.1) / (MII.total + 0.1))) %>% 
  dplyr::mutate(gene_biotype = factor(gene_biotype, levels = c("lnc1", "dormant", "histone", "protein_coding")), 
                gene_description = str_remove(gene_description, " \\[.*"),
                gene_description = replace(gene_description, is.na(gene_description), "")) %>% 
  dplyr::arrange(desc(gene_biotype))

# plot
cross_plot <-
  ggplot(plot_tb, aes(x = GV_ratio, y = MII_ratio, color = gene_biotype, size = gene_biotype, alpha = gene_biotype)) +
  geom_point(shape = 16) +
  scale_x_continuous(limits = c(-11, 11), breaks = seq(-10, 10, 2)) +
  scale_y_continuous(limits = c(-11, 11), breaks = seq(-10, 10, 2)) +
  scale_color_manual(values = c(lnc1 = "red3", dormant = "#1a75ff", histone = "black", protein_coding = "black")) +
  scale_size_manual(values = c(lnc1 = 4, dormant = 2.5, histone = 2.5, protein_coding = 2.5)) +
  scale_alpha_manual(values = c(lnc1 = 1, dormant = 1, histone = 0.2, protein_coding = 0.2)) + 
  guides(color = FALSE, size = FALSE, alpha = FALSE) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlab(str_c("log2 (GV polyA / GV total)")) +
  ylab(str_c("log2 (MII polyA / MII total)")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = 0.3),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c("ensembl.93.GRCm38.p6.20180919", "FPKM.polyA_total_ratio.GV_MII",
                                           str_remove(experiment_polyA_name, "_.*"),
                                           str_remove(experiment_total_name, "_.*"), "png", sep = ".")),
       plot = cross_plot,
       width = 10, height = 10)

# interactive plot
interactive_plot <-
  plotly::plot_ly(data = plot_tb,
                  x = ~GV_ratio,
                  y = ~MII_ratio,
                  text = ~paste("</br> GV polyA FPKM: ", round(GV.polyA, 3),
                                "</br> GV total FPKM: ", round(GV.total, 3),
                                "</br> MII polyA FPKM: ", round(MII.polyA, 3),
                                "</br> MII total FPKM: ", round(MII.total, 3),
                                "</br>", gene_id,
                                "</br>", gene_name,
                                "</br>", gene_description),
                  color = ~gene_biotype,
                  colors = c("red3", "#1a75ff", "gray60", "gray60"),
                  alpha = 0.75,
                  hoverinfo = "text") %>%
  add_markers() %>%
  layout(xaxis = list(title = "log2 (GV polyA / GV total)"),
         yaxis = list(title = "log2 (MII polyA / MII total)"))

# save as html widget
htmlwidgets::saveWidget(plotly::as_widget(interactive_plot),
                        file = file.path(outpath, str_c("ensembl.93.GRCm38.p6.20180919", "FPKM.polyA_total_ratio.GV_MII", 
                                                        str_remove(experiment_polyA_name, "_.*"), 
                                                        str_remove(experiment_total_name, "_.*"), "html", sep = ".")), 
                        selfcontained = T)
