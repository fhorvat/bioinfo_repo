### INFO: 
### DATE: Sun Oct 06 17:21:37 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_proteins_expression/results")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# annotation path
annotation_path <- file.path(inpath, "piRNA_proteins.gene_ids.with_notes.csv")

# coordinates path
coordinates_path <- file.path(inpath, "piRNA_proteins.coordinates.csv")

# datasets path
datasets_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_proteins_expression/datasets"

# FPKM table path
fpkm_paths <- list.files(datasets_path, ".*\\.FPKM_long\\.csv", full.names = T, recursive =  T)

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(annotation_path)

# read FPKM table
fpkm_tb_list <- purrr::map(fpkm_paths, function(path){
  
  # read FPKM
  fpkm_tb <- readr::read_csv(path)
  
  # get only GV for cow oocytes
  if(str_detect(path, "cow_oocyte"))fpkm_tb %<>% dplyr::filter(str_detect(sample_id, "GV"))
  
  # get only mature for cow testis
  if(str_detect(path, "cow_testis")) fpkm_tb %<>% dplyr::filter(str_detect(sample_id, "mature"))
  
  # get only WT for golden hamster
  if(str_detect(path, "golden_hamster_testis")) fpkm_tb %<>% dplyr::filter(str_detect(sample_id, "Mov10l_WT_adult"))

  # return 
  return(fpkm_tb)
  
}) %>% 
  set_names(., str_extract_all(fpkm_paths, "cow|mouse|rat|golden_hamster|human|testis|oocyte") %>% purrr::map(., str_c, collapse = ".") %>% unlist(.))

######################################################## MAIN CODE
### FPKM table
# create table
pirna_genes_tb <- 
  tibble(gene_name = c("Henmt1", "Mael", "Mov10l1", 
                       "Piwil1", "Piwil2", "Piwil3", "Piwil4", 
                       "Pld6", 
                       "Rnf17", "Spocd1", 
                       "Tdrd1", "Tdrd12", "Tdrd6", "Tdrd9", "Tdrkh", 
                       "Tex15"),
         other_name = c("HENMT1", "MAEL/Maelstorm", "MOV10L1", 
                        "MIWI", "MILI", "PIWIL3", "MIWI2", 
                        "PLD6/MITOPOLD/Zucchini", 
                        "TDRD4/RNF17", "SPOCD1", 
                        "TDRD1", "TDRD12", "TDRD6", "TDRD9", "TDRD2/TDRKH", 
                        "TEX15")) %>% 
  dplyr::arrange(gene_name)

# get FPKM of chosen genes
pirna_genes_fpkm <- purrr::map(names(fpkm_tb_list), function(animal_tissue){
  
  # get animal name
  animal <- str_remove(animal_tissue, "\\..*")
  
  # filter annotation
  annotation_filt <- 
    genes_info %>% 
    dplyr::select(gene_name, contains(animal)) %>% 
    set_names(c("gene_name", "gene_id")) %>% 
    dplyr::mutate(animal_tissue = animal_tissue)
  
  # filter fpkm table
  fpkm_filt <- 
    fpkm_tb_list[[animal_tissue]] %>% 
    dplyr::filter(gene_id %in% annotation_filt$gene_id) %>% 
    dplyr::select(gene_id, sample_id, fpkm) %>% 
    dplyr::group_by(gene_id) %>% 
    dplyr::summarize(mean = mean(fpkm), 
                     sd = sd(fpkm), 
                     length = length(fpkm)) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::mutate(SEM = sd / sqrt(length)) %>% 
    dplyr::right_join(., annotation_filt, by = "gene_id") %>% 
    dplyr::select(gene_name, gene_id, mean, sd, SEM, animal_tissue) %>% 
    dplyr::left_join(., pirna_genes_tb, by = "gene_name") %>% 
    dplyr::mutate(gene_name = factor(gene_name, levels = pirna_genes_tb$gene_name),
                  other_name = factor(other_name, levels = pirna_genes_tb$other_name))
  
  
  # return
  return(fpkm_filt)
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(animal_tissue = factor(animal_tissue, levels = str_c(c("cow", "human", "golden_hamster", "rat", "mouse"), 
                                                                     rep(c("testis", "oocyte"), each = 5), sep = ".")))

# transform to log10, split by oocyte and testis
pirna_plot_tb_list <- 
  pirna_genes_fpkm %>% 
  dplyr::mutate(mean = log10(mean + 0.001)) %>% 
  split(., str_extract(.$animal_tissue, "oocyte|testis"))

# get max and min values
pirna_plot_tb_limits <- 
  pirna_plot_tb_list %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::summarise(scale_max = max(mean, na.rm = T), 
                   scale_min = min(mean, na.rm = T))

# plot separately for oocyte and testis
purrr::map(names(pirna_plot_tb_list), function(tissue){
  
  # extract dataset
  pirna_plot_tb <- pirna_plot_tb_list[[tissue]]
  
  # heatmap ggplot2
  heat_plot <- 
    ggplot(pirna_plot_tb, aes(x = animal_tissue, y = gene_name)) + 
    geom_tile(aes(fill = mean), colour = "grey45") + 
    coord_equal() + 
    scale_fill_gradient2(low = "white", high = "black", na.value = "white", midpoint = -1, 
                         limits = c(pirna_plot_tb_limits$scale_min, pirna_plot_tb_limits$scale_max)) +
    theme(axis.text.x = element_text(size = 12, angle = 90, face = "bold", colour = "grey25", vjust = 0.5, hjust = 0), 
          axis.text.y = element_text(size = 12, face = "bold", colour = "grey25"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.position = "bottom", 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA, colour = NA), 
          axis.ticks = element_blank()) + 
    labs(x = "", 
         y = "", 
         fill = "") + 
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev(levels(pirna_genes_fpkm$gene_name))) 
  
  # save
  ggsave(filename = file.path(outpath, str_c("piRNA_pathway_genes.log_FPKM.animals", "with_human", tissue, "heatmap.png", sep = ".")), plot = heat_plot, width = 10, height = 10)
  
  # return
  return(tissue)
  
})


### save FPKM values
# long table to wide
fpkm_wide <- 
  pirna_genes_fpkm %>% 
  dplyr::select(gene_id, gene_name, other_name, animal_tissue, mean) %>% 
  tidyr::pivot_wider(id_cols = c("gene_name", "other_name"), names_from = animal_tissue, values_from = mean) %>% 
  dplyr::arrange(gene_name)

# save
readr::write_csv(fpkm_wide, path = file.path(outpath, "piRNA_proteins.FPKM.with_human.csv"))
readr::write_csv(pirna_genes_fpkm, path = file.path(outpath, "piRNA_proteins.FPKM.long.with_human.csv"))

# SEM
pirna_genes_fpkm %>% 
  dplyr::select(gene_id, gene_name, other_name, animal_tissue, SEM) %>% 
  tidyr::pivot_wider(id_cols = c("gene_name", "other_name"), names_from = animal_tissue, values_from = SEM) %>% 
  dplyr::arrange(gene_name) %T>% 
  readr::write_csv(., path = file.path(outpath, "piRNA_proteins.SEM.with_human.csv"))


# ### barplot
# # prepare data
# pirna_barplot_tb <-
#   pirna_genes_fpkm %>%
#   tidyr::separate(col = animal_tissue, into = c("animal", "tissue"), sep = "\\.")
# 
# # barplot ggplot2
# bar_plot <-
#   ggplot() +
#   geom_bar(data = pirna_barplot_tb,
#            mapping = aes(x = other_name, y = mean, fill = animal, group = interaction(animal, tissue)),
#            width = 0.8, stat = "identity",
#            position = position_dodge(width = 0.8)) +
#   geom_errorbar(data = pirna_barplot_tb,
#                 mapping = aes(x = other_name, ymin = mean - SEM, ymax = mean + SEM, color = animal, group = interaction(animal, tissue)),
#                 width = 0.4, size = 1,
#                 position = position_dodge(width = 0.8)) +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d() +
#   scale_x_discrete(drop = FALSE) +
#   theme_bw() +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank()) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "bottom")
# 
# # save
# ggsave(filename = file.path(outpath, "piRNA_pathway_genes.FPKM.animals.with_human.barplot.png"), plot = bar_plot, width = 15, height = 10)
# 
# 
