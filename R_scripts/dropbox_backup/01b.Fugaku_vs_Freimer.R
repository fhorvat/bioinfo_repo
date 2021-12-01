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

######################################################## READ DATA
# read FPKM tables
FPKM_total <- readr::read_csv(FPKM_total_path)
FPKM_polyA <- readr::read_csv(FPKM_polyA_path)

######################################################## MAIN CODE
# set lncRNA gene biotypes
lncrna_biotypes <- c("3prime_overlapping_ncRNA", "antisense", "macro_lncRNA", "sense_intronic", "sense_overlapping", "lincRNA", "bidirectional_promoter_lncRNA")

### tidy and join tables
# total table
FPKM_total_tidy <-
  FPKM_total %>% 
  dplyr::select(gene_id, GV.total = GV.WE, MII.total = MII.WE, gene_biotype, gene_name, gene_description)

# polyA table
FPKM_polyA_tidy <-
  FPKM_polyA %>% 
  dplyr::select(gene_id, GV.polyA = GV, MII.polyA = MII, gene_biotype)

# join tables
FPKM_tb <- 
  left_join(FPKM_total_tidy, FPKM_polyA_tidy, by = c("gene_id", "gene_biotype")) %>% 
  dplyr::filter(gene_biotype %in% c(lncrna_biotypes, "protein_coding")) %>% 
  dplyr::mutate(gene_biotype = ifelse(gene_biotype %in% lncrna_biotypes, "lncRNA", gene_biotype), 
                gene_biotype = replace(gene_biotype, gene_id == "ENSMUSG00000110001", "lnc1")) %>% 
  dplyr::select(gene_id, GV.polyA, GV.total, MII.polyA, MII.total, gene_biotype, gene_name, gene_description)


### construct and save cross-hair plot
# table for plot
plot_tb <- 
  FPKM_tb %>% 
  dplyr::mutate(GV_ratio = log2((GV.polyA + 0.1) / (GV.total + 0.1)), 
                MII_ratio = log2((MII.polyA + 0.1) / (MII.total + 0.1))) %>% 
  dplyr::mutate(gene_biotype = factor(gene_biotype, levels = c("lnc1", "lncRNA", "protein_coding")), 
                gene_description = str_remove(gene_description, " \\[.*"),
                gene_description = replace(gene_description, is.na(gene_description), "")) %>% 
  dplyr::arrange(desc(gene_biotype))

# plot
cross_plot <-
  ggplot(plot_tb, aes(x = GV_ratio, y = MII_ratio, color = gene_biotype, size = gene_biotype)) +
  geom_point() +
  scale_x_continuous(limits = c(-10.5, 10.5), breaks = seq(-10, 10, 2)) +
  scale_y_continuous(limits = c(-10.5, 10.5), breaks = seq(-10, 10, 2)) +
  scale_color_manual(values = c(protein_coding = "black", lncRNA = "#1a75ff", lnc1 = "red3")) +
  scale_size_manual(values = c(protein_coding = 2.5, lncRNA = 2.5, lnc1 = 4)) +
  guides(color = FALSE, size = FALSE) +
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
                  colors = c("red3", "#1a75ff", "black"),
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
