### INFO: produce scatterplots of ratio of expression in Fugaku and CNOT6L data
### DATE: Sun Mar 11 00:36:34 2018
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

library(VennDiagram)
library(geneplotter)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "mutate_cond.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set outpath
outpath <- getwd()

# ENSEMBL annotated genes info path
ensembl_genes_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.89.GRCm38.p5.20180615.UCSCseqnames.geneInfo.csv"

######################################################## READ DATA
# read info about all ENSEMBL annotated genes
ensembl_genes_info <- readr::read_csv(ensembl_genes_path)

# CNOT6L FPKM
fpkm_CNOT6L <- readr::read_csv(file = file.path(outpath, "ensembl.GRCm38.89.CNOT6L.avgFPKM.csv"))

# Fugaku FPKM
fpkm_Fugaku <- readr::read_csv(file = file.path(outpath, "ensembl.GRCm38.89.Fugaku.FPKM.csv"))

# Yu FPKM
fpkm_Yu <- 
  readr::read_csv(file = file.path(outpath, "Mus_musculus.GRCm38.89.20180305.Yu2016.avgFPKM.csv")) %>% 
  data.table::setnames(., old = 2:ncol(.), new = str_c(colnames(.)[2:ncol(.)], ".Yu") %>% stringr::str_remove(., pattern = ".PE"))

# CNOT6L significantly diff. exp. genes
results_CNOT6L <- 
  list.files(file.path(outpath, "results"), "diffExp.CNOT6L.*.signif.csv", full.names = T) %>% 
  lapply(., readr::read_csv) %>% 
  dplyr::bind_rows(.)

# Ma 2013 significantly diff. exp. genes
results_Ma <- readr::read_csv(file.path(outpath, "affy.mouse4302.Ma2013.controlVsDcp1a_Dcp2.upregulated.csv"))

# Yu 2016 significantly upregulated genes
results_Yu <- 
  list.files(file.path(outpath, "results"), "diffExp.Yu2016.*.signif.csv", full.names = T) %>% 
  lapply(., readr::read_csv) %>% 
  dplyr::bind_rows(.)

######################################################## MAIN CODE
# get gene_id of protein coding genes
protein_coding <- 
  ensembl_genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# filter results table
results <- 
  results_CNOT6L %>% 
  dplyr::mutate(stage = str_extract(comparison, "1C|MII|GV"), 
                regulation = ifelse(log2FoldChange > 0, "up", "down")) %>% 
  dplyr::select(gene_id, stage, regulation, logFC_CNOT6L = log2FoldChange, gene_description, 
                MII_KO_FPKM = MII_KO, MII_WT_FPKM = MII_WT)

# filter Yu results table
results_Yu %<>% 
  dplyr::mutate(stage = str_extract(comparison, "1C|MII|GV"), 
                regulation = ifelse(log2FoldChange > 0, "up", "down")) %>% 
  dplyr::select(gene_id, stage, regulation, logFC_Yu = log2FoldChange, gene_description) %>% 
  dplyr::filter(gene_id %in% protein_coding) 

# join tables, filter only protein coding
fpkm_df <- 
  dplyr::left_join(fpkm_CNOT6L, fpkm_Fugaku, by = "gene_id") %>% 
  dplyr::left_join(., fpkm_Yu, by = "gene_id") %>% 
  dplyr::filter(gene_id %in% protein_coding) %>% 
  data.table::setnames(., 2:ncol(.), str_c("s.", colnames(.)[2:ncol(.)]))

# filter Ma
results_Ma %<>%
  dplyr::filter(gene_id %in% protein_coding) %>% 
  dplyr::rename(logFC_Ma = logFC)


###### plot CNOT6L/Fugaku ratio scatter plots
### column names
## CNOT6L: s.1C_KO, s.1C_WT, s.GV_KO, s.GV_WT, s.MII_KO, s.MII_WT

## Fugaku:
# s.GV.WE s.MII.WE s.1C.WE s.2C.WE s.4C.WE s.Blast.WE s.Molura.WE
# s.1C.PA s.MII.PA 

### plot ratio scatter plots
## plot 1 
# - x: s.MII_WT / s.GV_WT  (CNOT6L)
# - y: s.MII.WE / s.GV.WE (Fugaku)
# - signif: MII (CNOT6L)

## plot 2
# - x: s.1C_WT / s.MII_WT (CNOT6L)
# - y: s.1C.WE / s.MII.WE (Fugaku)
# - signif: 1C (CNOT6L)

## plot 3:
# - x: s.1C.PA / s.MII.PA (Fugaku)
# - y: s.1C.WE / s.MII.WE (Fugaku)
# - signif: 1C (CNOT6L)

### plot
# set variables
varible_list <- list(x1 = c("s.MII_WT", "s.1C_WT", "s.1C.PA"),
                     x2 = c("s.GV_WT", "s.MII_WT", "s.MII.PA"),
                     y1 = c("s.MII.WE", "s.1C.WE", "s.1C.WE"),
                     y2 = c("s.GV.WE", "s.MII.WE", "s.MII.WE"),
                     x_set = c("CNOT6L", "CNOT6L", "Fugaku"),
                     y_set = c("Fugaku", "Fugaku", "Fugaku"),
                     signif_stage = c("MII", "1C", "1C"))

# plot in loop
for(n in 1:length(varible_list[[1]])){
  
  # subset variables
  x1 <- varible_list$x1[n]
  x2 <- varible_list$x2[n]
  y1 <- varible_list$y1[n]
  y2 <- varible_list$y2[n]
  x_set <- varible_list$x_set[n]
  y_set <- varible_list$y_set[n]
  signif_stage <- varible_list$signif_stage[n]
  
  # create ratio data.frame
  ratio_df <-
    fpkm_df %>% 
    dplyr::mutate_(x1 = x1, x2 = x2, y1 = y1, y2 = y2) %>%
    dplyr::select(x1, x2, y1, y2, gene_id) %>%
    dplyr::filter_at(.vars = vars(x2, y2), .vars_predicate = all_vars(. > 1)) %>%
    dplyr::mutate(x = log2(x1 / x2),
                  y = log2(y1 / y2)) %>%
    dplyr::select(gene_id, x, y) %>%
    dplyr::filter(complete.cases(.),
                  !is.infinite(x),
                  !is.infinite(y))
  
  # plot
  plot_ratio <-
    ratio_df %>%
    dplyr::left_join(., results %>% dplyr::filter(stage == signif_stage) %>% dplyr::select(-stage), by = "gene_id") %>%
    dplyr::mutate(regulation = replace(regulation, is.na(regulation), "no"), 
                  regulation = factor(regulation, levels = c("no", "down", "up"))) %>%
    dplyr::arrange(regulation) %>%
    ggplot(data = ., aes(x = x, y = y, color = regulation)) +
    geom_point(size = 2.5) +
    scale_x_continuous(limits = c(- 8.1, 8.1), breaks = seq(-8, 8, 2)) +
    scale_y_continuous(limits = c(- 8.1, 8.1), breaks = seq(-8, 8, 2)) +
    scale_color_manual(values = c(up = "red3", down = "gray60", no = "gray60"), 
                       breaks = c("up", "no")) +
    guides(color = FALSE) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    xlab(str_c("log2 (", str_replace(x1, "s.", ""), " / ", str_replace(x2, "s.", ""), ") ", x_set)) +
    ylab(str_c("log2 (", str_replace(y1, "s.", ""), " / ", str_replace(y2, "s.", ""), ") ", y_set)) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
          axis.title.y = element_text(size = 15, vjust = 0.3),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save plot
  ggsave(filename = file.path(outpath, "results", str_c("scatterRatio.comp.",
                                                        str_replace(x1, "s.", ""), "_", str_replace(x2, "s.", ""), ".", x_set, ".vs.",
                                                        str_replace(y1, "s.", ""), "_", str_replace(y2, "s.", ""), ".", y_set, 
                                                        ".filterEarlyStage.1FPKM.png")), 
         plot = plot_ratio, width = 7.5, height = 15)
  
}
###### 


###### comparison with Dcp KOs from Ma 2013 paper
# set variables
varible_list <- list(x1 = c("s.MII_WT", "s.1C_WT", "s.1C.PA"),
                     x2 = c("s.GV_WT", "s.MII_WT", "s.MII.PA"),
                     y1 = c("s.MII.WE", "s.1C.WE", "s.1C.WE"),
                     y2 = c("s.GV.WE", "s.MII.WE", "s.MII.WE"),
                     x_set = c("CNOT6L", "CNOT6L", "Fugaku"),
                     y_set = c("Fugaku", "Fugaku", "Fugaku"),
                     signif_stage = c("MII", "MII", "MII"))


# plot in loop
for(n in 1:length(varible_list[[1]])){
  
  # subset variables
  x1 <- varible_list$x1[n]
  x2 <- varible_list$x2[n]
  y1 <- varible_list$y1[n]
  y2 <- varible_list$y2[n]
  x_set <- varible_list$x_set[n]
  y_set <- varible_list$y_set[n]
  signif_stage <- varible_list$signif_stage[n]
  
  # create ratio data.frame
  ratio_df <-
    fpkm_df %>% 
    dplyr::mutate_(x1 = x1, x2 = x2, y1 = y1, y2 = y2) %>%
    dplyr::select(x1, x2, y1, y2, gene_id) %>%
    dplyr::filter_at(.vars = vars(x2, y2), .vars_predicate = all_vars(. > 1)) %>%
    dplyr::mutate(x = log2(x1 / x2),
                  y = log2(y1 / y2)) %>%
    dplyr::select(gene_id, x, y) %>%
    dplyr::filter(complete.cases(.),
                  !is.infinite(x),
                  !is.infinite(y)) %>%
    dplyr::left_join(., results %>% dplyr::filter(stage == signif_stage) %>% dplyr::select(-stage), by = "gene_id") %>%
    dplyr::mutate(regulation = replace(regulation, is.na(regulation), "no"),
                  Ma_upregulated = (gene_id %in% results_Ma$gene_id),
                  CNOT6L_upregulated = (regulation == "up")) %>%
    mutate_cond(Ma_upregulated, upregulation = "Ma") %>%
    mutate_cond(CNOT6L_upregulated, upregulation = "CNOT6L") %>%
    mutate_cond(Ma_upregulated & CNOT6L_upregulated, upregulation = "both") %>%
    mutate_cond(!(Ma_upregulated | CNOT6L_upregulated), upregulation = "neither") %>%
    dplyr::mutate(upregulation = factor(upregulation, levels = c("neither", "CNOT6L", "Ma", "both"))) %>%
    dplyr::arrange(upregulation)
  
  # plot
  plot_ratio <-
    ggplot(data = ratio_df, aes(x = x, y = y, color = upregulation)) +
    geom_point(size = 1.25) +
    scale_x_continuous(limits = c(- 8.1, 8.1), breaks = seq(-8, 8, 2)) +
    scale_y_continuous(limits = c(- 8.1, 8.1), breaks = seq(-8, 8, 2)) +
    scale_color_manual(values = c(neither = "gray60", CNOT6L = "red3", Ma = "black", both = "#1a75ff")) +
    scale_alpha_manual(values = c(neither = 0.8, CNOT6L = 0.3, Ma = 0.8, both = 1)) +
    guides(color = FALSE) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    xlab(str_c("log2 (", str_replace(x1, "s.", ""), " / ", str_replace(x2, "s.", ""), ") ", x_set)) +
    ylab(str_c("log2 (", str_replace(y1, "s.", ""), " / ", str_replace(y2, "s.", ""), ") ", y_set)) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
          axis.title.y = element_text(size = 15, vjust = 0.3),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save plot
  ggsave(filename = file.path(outpath, "results", str_c("scatterRatio.comp.",
                                                        str_replace(x1, "s.", ""), "_", str_replace(x2, "s.", ""), ".", x_set, ".vs.",
                                                        str_replace(y1, "s.", ""), "_", str_replace(y2, "s.", ""), ".", y_set, 
                                                        ".Ma.filterEarlyStage.1FPKM.png")), 
         plot = plot_ratio, width = 7.5, height = 15)
}
###### 


###### plot CNOT6L ratio scatter plots
# on Y will be WT MII/GV and on X MII-KO/MII-WT
# on Y will be WT 1C/MII and on X 1C-KO/1C-WT
# and the blue and red genes the up and down in MII and 1C, respectively

# plot 1
# - x: s.MII_KO / s.MII_WT (CNOT6L)
# - y: s.MII_WT / s.GV_WT (CNOT6L)
# - signif: MII (CNOT6L)

# plot 2
# - x: s.1C_KO / s.1C_WT (CNOT6L)
# - y: s.1C_WT / s.MII_WT (CNOT6L)
# - signif: 1C (CNOT6L)

### plot
# set variables
varible_list <- list(x1 = c("s.MII_KO", "s.1C_KO"),
                     x2 = c("s.MII_WT", "s.1C_WT"),
                     y1 = c("s.MII_WT", "s.1C_WT"),
                     y2 = c("s.GV_WT", "s.MII_WT"),
                     x_set = c("CNOT6L", "CNOT6L"),
                     y_set = c("CNOT6L", "CNOT6L"),
                     signif_stage = c("MII", "1C"))

# plot in loop
for(n in 1:length(varible_list[[1]])){
  
  # subset variables
  x1 <- varible_list$x1[n]
  x2 <- varible_list$x2[n]
  y1 <- varible_list$y1[n]
  y2 <- varible_list$y2[n]
  x_set <- varible_list$x_set[n]
  y_set <- varible_list$y_set[n]
  signif_stage <- varible_list$signif_stage[n]
  
  # create ratio data.frame
  ratio_df <-
    fpkm_df %>% 
    dplyr::mutate_(x1 = x1, x2 = x2, y1 = y1, y2 = y2) %>%
    dplyr::select(x1, x2, y1, y2, gene_id) %>% 
    dplyr::filter_at(.vars = vars(x2, y2), .vars_predicate = all_vars(. > 1)) %>%
    dplyr::mutate(x = log2(x1 / x2),
                  y = log2(y1 / y2)) %>%
    dplyr::select(gene_id, x, y) %>%
    dplyr::filter(complete.cases(.),
                  !is.infinite(x),
                  !is.infinite(y))
  
  # create plot data.frame, plot
  plot_ratio <-
    ratio_df %>%
    dplyr::left_join(., results %>% dplyr::filter(stage == signif_stage) %>% dplyr::select(-stage), by = "gene_id") %>%
    dplyr::mutate(regulation = replace(regulation, is.na(regulation), "no"), 
                  regulation = factor(regulation, levels = c("no", "down", "up"))) %>%
    dplyr::arrange(regulation) %>%
    ggplot(data = ., aes(x = x, y = y, color = regulation)) +
    geom_point(size = 5) +
    scale_x_continuous(limits = c(-8.1, 8.1), breaks = seq(-8, 8, 2)) +
    scale_y_continuous(limits = c(-8.1, 8.1), breaks = seq(-8, 8, 2)) +
    scale_color_manual(values = c(up = "red3", down = "#1a75ff", no = "gray60")) +
    guides(color = FALSE) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    xlab(str_c("log2 (", str_replace(x1, "s.", ""), " / ", str_replace(x2, "s.", ""), ") ", x_set)) +
    ylab(str_c("log2 (", str_replace(y1, "s.", ""), " / ", str_replace(y2, "s.", ""), ") ", y_set)) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
          axis.title.y = element_text(size = 15, vjust = 0.3),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save plot
  ggsave(filename = file.path(outpath, "results", str_c("scatterRatio.CNOT6L.",
                                                        str_replace(x1, "s.", ""), "_", str_replace(x2, "s.", ""), ".vs.",
                                                        str_replace(y1, "s.", ""), "_", str_replace(y2, "s.", ""),
                                                        ".filterWT.1FPKM.png")), 
         plot = plot_ratio, width = 15, height = 15)
  
}
###### 


###### plot BTG4 ratio scatter plots
# on Y will be Yu WT MII/GV and on X Yu MII-KO/MII-WT
# on Y will be CNOT6L WT MII/GV and on X Yu MII-KO/MII-WT
# and the blue and red genes the up and down in MII and MII, respectively

# plot 1
# - x: s.MII_KO.Yu / s.MII_WT.Yu (Yu)
# - y: s.MII_WT.Yu / s.GV_WT.Yu (Yu)
# - signif: MII (Yu)

# plot 2
# - x: s.MII_KO.Yu / s.MII_WT.Yu (Yu)
# - y: s.MII_WT / s.GV_WT (CNOT6L)
# - signif: MII (Yu)

### plot
# set variables
varible_list <- list(x1 = c("s.MII_KO.Yu", "s.MII_KO.Yu"),
                     x2 = c("s.MII_WT.Yu", "s.MII_WT.Yu"),
                     y1 = c("s.MII_WT.Yu", "s.MII_WT"),
                     y2 = c("s.GV_WT.Yu", "s.GV_WT"),
                     x_set = c("Yu", "Yu"),
                     y_set = c("Yu", "CNOT6L"),
                     signif_stage = c("MII", "MII"))

# plot in loop
for(n in 1:length(varible_list[[1]])){
  
  # subset variables
  x1 <- varible_list$x1[n]
  x2 <- varible_list$x2[n]
  y1 <- varible_list$y1[n]
  y2 <- varible_list$y2[n]
  x_set <- varible_list$x_set[n]
  y_set <- varible_list$y_set[n]
  signif_stage <- varible_list$signif_stage[n]
  
  # create ratio data.frame
  ratio_df <-
    fpkm_df %>% 
    dplyr::mutate_(x1 = x1, x2 = x2, y1 = y1, y2 = y2) %>%
    dplyr::select(x1, x2, y1, y2, gene_id) %>% 
    dplyr::filter_at(.vars = vars(x2, y2), .vars_predicate = all_vars(. > 1)) %>%
    dplyr::mutate(x = log2(x1 / x2),
                  y = log2(y1 / y2)) %>%
    dplyr::select(gene_id, x, y) %>%
    dplyr::filter(complete.cases(.),
                  !is.infinite(x),
                  !is.infinite(y)) %>% 
    dplyr::left_join(., results_Yu %>% dplyr::filter(stage == signif_stage) %>% dplyr::select(-stage), by = "gene_id") %>%
    dplyr::mutate(regulation = replace(regulation, is.na(regulation), "no")) %>% 
    dplyr::mutate(regulation = factor(regulation, levels = c("no", "down", "up"))) %>%
    dplyr::arrange(regulation)
  
  # x <- 
  #   ratio_df %>% 
  #   dplyr::filter(regulation == "up", x < 0) %$%
  #   gene_id

  # create plot data.frame, plot
  plot_ratio <-
    ggplot(data = ratio_df, aes(x = x, y = y, color = regulation)) +
    geom_point(size = 1) +
    scale_x_continuous(limits = c(-8.1, 8.1), breaks = seq(-8, 8, 2)) +
    scale_y_continuous(limits = c(-8.1, 8.1), breaks = seq(-8, 8, 2)) +
    scale_color_manual(values = c(up = "red3", down = "#1a75ff", no = "gray60")) +
    guides(color = FALSE) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    xlab(str_c("log2 (", str_replace_all(x1, "s.|.Yu", ""), " / ", str_replace_all(x2, "s.|.Yu", ""), ") ", x_set)) +
    ylab(str_c("log2 (", str_replace_all(y1, "s.|.Yu", ""), " / ", str_replace_all(y2, "s.|.Yu", ""), ") ", y_set)) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
          axis.title.y = element_text(size = 15, vjust = 0.3),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save plot
  ggsave(filename = file.path(outpath, "results", str_c("scatterRatio.BTG4.",
                                                        str_replace_all(x1, "s.|.Yu", ""), "_", str_replace_all(x2, "s.|.Yu", ""), ".", x_set, ".vs.",
                                                        str_replace_all(y1, "s.|.Yu", ""), "_", str_replace_all(y2, "s.|.Yu", ""), ".", y_set, 
                                                        ".filterEarlyStage.1FPKM.png")), 
         plot = plot_ratio, width = 15, height = 15)
  
}
###### 


# ###### comparison with Dcp KOs from Ma 2013 paper
# x1 <- "s.MII_WT"
# x2 <- "s.GV_WT"
# y1 <- "s.MII.WE"
# y2 <- "s.GV.WE"
# x_set <- "CNOT6L"
# y_set <- "Fugaku"
# signif_stage <- "MII"
# 
# # create ratio data.frame
# ratio_df <-
#   fpkm_df %>%
#   dplyr::mutate_(x1 = x1, x2 = x2, y1 = y1, y2 = y2) %>%
#   dplyr::select(x1, x2, y1, y2, gene_id) %>%
#   # dplyr::filter_at(.vars = vars(y1, y2), .vars_predicate = any_vars(. > 1)) %>%
#   # dplyr::filter_at(.vars = vars(x1, x2), .vars_predicate = any_vars(. > 1)) %>%
#   dplyr::mutate(x = log2(x1 / x2),
#                 y = log2(y1 / y2)) %>%
#   dplyr::select(gene_id, x, y) %>%
#   dplyr::filter(complete.cases(.),
#                 !is.infinite(x),
#                 !is.infinite(y))
# 
# # create plot data.frame
# plot_ratio_df <-
#   ratio_df %>%
#   dplyr::left_join(., results %>% dplyr::filter(stage == signif_stage) %>% dplyr::select(-stage), by = "gene_id") %>%
#   dplyr::mutate(regulation = replace(regulation, is.na(regulation), "no"),
#                 Ma_upregulated = (gene_id %in% results_Ma$gene_id),
#                 CNOT6L_upregulated = (regulation == "up")) %>%
#   mutate_cond(Ma_upregulated, upregulation = "Ma") %>%
#   mutate_cond(CNOT6L_upregulated, upregulation = "CNOT6L") %>%
#   mutate_cond(Ma_upregulated & CNOT6L_upregulated, upregulation = "both") %>%
#   mutate_cond(!(Ma_upregulated | CNOT6L_upregulated), upregulation = "neither") %>%
#   dplyr::mutate(upregulation = factor(upregulation, levels = c("neither", "CNOT6L", "Ma", "both"))) %>%
#   dplyr::arrange(upregulation)
# 
#
# ### get list of genes for which both KO have up-regulating effect
# plot_ratio_both <- 
#   plot_ratio_df %>% 
#   dplyr::filter(upregulation == "both") %>% 
#   dplyr::select(gene_id) %>% 
#   dplyr::left_join(., results_Ma, by = "gene_id") %>% 
#   dplyr::left_join(., results %>% dplyr::filter(stage == signif_stage), by = "gene_id") %>% 
#   dplyr::mutate(gene_description = str_replace(gene_description, " \\[Source.*", "")) %>% 
#   dplyr::group_by(gene_symbol) %>% 
#   dplyr::filter(logFC_Ma == max(logFC_Ma)) %>% 
#   dplyr::ungroup(.) %>%
#   dplyr::select(gene_id, gene_symbol, logFC_Ma, logFC_CNOT6L, MII_WT_FPKM, MII_KO_FPKM, gene_description) %T>% 
#   readr::write_csv(., file.path(outpath, "results", "Ma_vs_CNOT6L.significant.upregulated.both.csv"))
# 
# ###### 

# ###### control
# ### column names
# ## CNOT6L: 
# # s.GV_KO s.GV_WT s.MII_KO s.MII_WT s.1C_KO s.1C_WT
# 
# ## Fugaku:
# # s.GV.WE s.MII.WE s.1C.WE s.2C.WE s.4C.WE s.Blast.WE s.Molura.WE
# # s.1C.PA s.MII.PA 
# 
# # set variables
# varible_list <- list(x = c("s.MII.WE", "s.1C.WE", "s.GV_WT", "s.MII_WT", "s.1C_WT", "s.MII_WT", "s.1C_WT"),
#                      y = c("s.MII.PA", "s.1C.PA", "s.GV.WE", "s.MII.WE", "s.1C.WE", "s.MII.PA", "s.1C.PA"),
#                      x_set = c("Fugaku", "Fugaku", "CNOT6L", "CNOT6L", "CNOT6L", "CNOT6L", "CNOT6L"),
#                      y_set = c("Fugaku", "Fugaku", "Fugaku", "Fugaku", "Fugaku", "Fugaku", "Fugaku"),
#                      signif_stage = c("MII", "1C", "GV", "MII", "1C", "MII", "1C"))
# 
# # plot in loop
# for(n in 1:length(varible_list[[1]])){
#   
#   # subset variables
#   x <- varible_list$x[n]
#   y <- varible_list$y[n]
#   x_set <- varible_list$x_set[n]
#   y_set <- varible_list$y_set[n]
#   signif_stage <- varible_list$signif_stage[n]
#   
#   # create ratio data.frame
#   ratio_df <-
#     fpkm_df %>% 
#     dplyr::mutate_(x = x, y = y) %>%
#     dplyr::select(gene_id, x, y) %>% 
#     # dplyr::filter_at(.vars = vars(y, y), .vars_predicate = any_vars(. > 1)) %>%
#     dplyr::mutate(x = log2(x + 1), 
#                   y = log2(y + 1)) %>% 
#     dplyr::filter(complete.cases(.),
#                   !is.infinite(x),
#                   !is.infinite(y))
#   
#   # create plot data.frame, plot
#   plot_ratio <-
#     ratio_df %>%
#     dplyr::left_join(., results %>% dplyr::filter(stage == signif_stage) %>% dplyr::select(-stage), by = "gene_id") %>%
#     dplyr::mutate(regulation = replace(regulation, is.na(regulation), "no"), 
#                   regulation = factor(regulation, levels = c("no", "down", "up"))) %>%
#     dplyr::arrange(regulation) %>%
#     ggplot(data = ., aes(x = x, y = y, color = regulation)) +
#     geom_point(size = 5) +
#     # scale_x_continuous(limits = c(-8.1, 8.1), breaks = seq(-8, 8, 2)) +
#     # scale_y_continuous(limits = c(-8.1, 8.1), breaks = seq(-8, 8, 2)) +
#     scale_color_manual(values = c(up = "red3", down = "#1a75ff", no = "gray60")) +
#     guides(color = FALSE) +
#     geom_vline(xintercept = 0) +
#     geom_hline(yintercept = 0) +
#     xlab(str_c("log2 (", str_replace(x, "s.", ""), ") ", x_set)) +
#     ylab(str_c("log2 (", str_replace(y, "s.", ""), ") ", y_set)) +
#     theme_bw() +
#     theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
#           axis.title.y = element_text(size = 15, vjust = 0.3),
#           axis.text.x = element_text(size = 15),
#           axis.text.y = element_text(size = 15),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank())
#   
#   # save plot
#   ggsave(filename = file.path(outpath, "results", str_c("scatterRatio.control.",
#                                                         str_replace(x, "s.", ""), ".", x_set, ".vs.",
#                                                         str_replace(y, "s.", ""), ".", y_set, ".",
#                                                         ".filterNone.png")), 
#          plot = plot_ratio, width = 15, height = 15)
# }
# ######