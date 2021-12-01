### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/comparison_with_mouse/Gan_2013_NatCommun_GSE35005")

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

library(scales)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames\\.geneInfo\\.csv$"), full.names = T)

# mouse expression path
mouse_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Gan_2013_NatCommun_GSE35005/Analysis/expression"
mouse_tb_path <- list.files(mouse_path, str_c("ensembl\\.", ensembl_version, ".*\\.FPKM_mean\\.csv.*$"), full.names = T)

# hamster expression path
hamster_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.RNAseq/Analysis/expression"
hamster_tb_path <- list.files(hamster_path, str_c("ensembl\\.", ensembl_version, ".*\\.FPKM_mean\\.csv.*$"), full.names = T)

# hamster-mouse orthologs path
ortho_path <- file.path(inpath, "..", "ensembl.99.mouse_vs_goldHamster.one2one_homologs.csv")

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read mouse expression
mouse_tb <- readr::read_csv(mouse_tb_path)

# read hamster expression
hamster_tb <- readr::read_csv(hamster_tb_path)

# read orthologs
ortho_tb <- readr::read_csv(ortho_path)

######################################################## MAIN CODE
### filter and join expression tables based on orthology
# mouse
mouse_ortho_fpkm <- 
  mouse_tb %>% 
  dplyr::left_join(., ortho_tb, by = c("gene_id" = "mouse_gene_id")) %>% 
  dplyr::select(mouse_gene_name, SG_A:round_ST) %>% 
  dplyr::filter(!is.na(mouse_gene_name)) %>% 
  magrittr::set_colnames(., c("mouse_gene_name", str_c("mouse.", colnames(.)[2:ncol(.)])))

# hamster
hamster_ortho_fpkm <- 
  hamster_tb %>% 
  dplyr::left_join(., ortho_tb, by = c("gene_id" = "hamster_gene_id")) %>% 
  dplyr::select(mouse_gene_name, Mov10l_WT_13dpp:Mov10l_WT_adult, Mov10l_KO_13dpp:Mov10l_WT_adult) %>% 
  dplyr::filter(!is.na(mouse_gene_name)) %>% 
  magrittr::set_colnames(., str_remove(colnames(.), "Mov10l_")) %>% 
  magrittr::set_colnames(., c("mouse_gene_name", str_c("hamster.", colnames(.)[2:ncol(.)], "_testis")))

# joined
ortho_fpkm <- 
  left_join(mouse_ortho_fpkm, hamster_ortho_fpkm, by = "mouse_gene_name") %>% 
  dplyr::rename(gene_name = "mouse_gene_name")

### calculate and plot correlation 
# get correlation
ortho_cc <- 
  ortho_fpkm %>% 
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "gene_name") %>% 
  as.matrix(.) %>% 
  cor(., method = "spearman") %>% 
  as_tibble(., rownames = "sample_mouse") %>% 
  tidyr::pivot_longer(-sample_mouse, names_to = "sample_hamster", values_to = "correlation") %>% 
  dplyr::filter(sample_mouse != "mouse.Sertoli") %>% 
  dplyr::filter(str_detect(sample_mouse, "mouse"), 
                str_detect(sample_hamster, "hamster")) %>% 
  dplyr::mutate(sample_mouse = factor(sample_mouse, levels = str_c("mouse.", c("primitive_SG_A", "SG_A", "SG_B",
                                                                                   "leptotene_SC", "pachytene_SC",
                                                                                   "round_ST", "elongative_ST"))), 
                sample_hamster = factor(sample_hamster, levels = str_c("hamster.", c("WT_13dpp", "WT_21dpp", "WT_adult",
                                                                                     "KO_13dpp", "KO_21dpp", "KO_adult"), "_testis")))
# plot using ggplot2
heat_plot <- 
  ggplot(ortho_cc, aes(x = sample_hamster, y = sample_mouse)) + 
  geom_tile(aes(fill = correlation), colour = "grey45") + 
  coord_equal() + 
  scale_fill_gradient(low = "white", high = "black") +
  # scale_fill_viridis_c() +
  theme(axis.text.x = element_text(size = 12, angle = 90, face = "bold", colour = "grey25", vjust = 0.5, hjust = 0), 
        axis.text.y = element_text(size = 12, face = "bold", colour = "grey25"), 
        legend.title = element_text(size = 10, face = "bold"), 
        legend.position = "bottom", 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = NA), 
        axis.ticks = element_blank()) + 
  labs(x = "", 
       y = "", 
       fill = "Pearsons's Correlation") + 
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(ortho_cc$sample_mouse))) 

# save
ggsave(plot = heat_plot, 
       filename = file.path(outpath, str_c("ensembl", ensembl_version, "hamster_vs_mouse.all_genes.FPKM.correlation.png", sep = ".")), 
       width = 10, height = 10)