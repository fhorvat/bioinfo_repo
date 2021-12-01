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
genes_info <- 
  readr::read_csv(annotation_path) %>% 
  dplyr::mutate(gene_id.bull_testis = gene_id.cow)

# read FPKM table
fpkm_tb_list <- 
  purrr::map(fpkm_paths, function(path){
    
    # read FPKM
    fpkm_tb <- readr::read_csv(path)
    
    # get only GV for cow
    if(str_detect(path, "cow")){
      
      fpkm_tb %<>% 
        dplyr::filter(str_detect(sample_id, "GV"))
      
    }
    
    # get only 
    if(str_detect(path, "bull")){
      
      fpkm_tb %<>% 
        dplyr::filter(str_detect(sample_id, "mature"))
      
    }
    
    # return 
    return(fpkm_tb)
    
  }) %>% 
  set_names(., str_extract(fpkm_paths, "cow|mouse|rat|golden_hamster|bull_testis"))

######################################################## MAIN CODE
### FPKM table
# create table
pirna_genes_tb <- 
  tibble(gene_name = c("Piwil2", "Piwil1", "Piwil4", "Piwil3",
                       "Mov10l1", "Ddx4", "Mael", 
                       "Pld6", "Henmt1",
                       "Tdrd1", "Tdrkh", "Rnf17", "Tdrd6", "Tdrd9", "Tdrd12", 
                       "Prmt5", "Wdr77", "Kif17"),
         other_name = c("MILI", "MIWI", "MIWI2", "PIWIL3",
                        "MOV10L1", "DDX4/VASE", "MAEL/Maelstorm", 
                        "PLD6/MITOPOLD/Zucchini", "HENMT1", 
                        "TDRD1", "TDRD2/TDRKH", "TDRD4/RNF17", "TDRD6", "TDRD9", "TDRD12", 
                        "PRMT5/Dart5/Capsuleen", "WDR77", "KIF17"))

# get FPKM of chosen genes
pirna_genes_fpkm <- 
  purrr::map(names(fpkm_tb_list), function(animal){
    
    # filter annotation
    annotation_filt <- 
      genes_info %>% 
      dplyr::select(gene_name, contains(animal)) %>% 
      set_names(c("gene_name", "gene_id")) %>% 
      dplyr::mutate(animal = animal)
    
    # filter fpkm table
    fpkm_filt <- 
      fpkm_tb_list[[animal]] %>% 
      dplyr::filter(gene_id %in% annotation_filt$gene_id) %>% 
      dplyr::select(gene_id, sample_id, fpkm) %>% 
      group_by(gene_id) %>% 
      summarize(mean = mean(fpkm), 
                sd = sd(fpkm), 
                length = length(fpkm)) %>% 
      dplyr::ungroup(.) %>% 
      dplyr::mutate(SEM = sd / sqrt(length)) %>% 
      dplyr::right_join(., annotation_filt, by = "gene_id") %>% 
      dplyr::select(gene_name, gene_id, mean, sd, SEM, animal) %>% 
      dplyr::left_join(., pirna_genes_tb, by = "gene_name") %>% 
      dplyr::mutate(gene_name = factor(gene_name, levels = pirna_genes_tb$gene_name),
                    other_name = factor(other_name, levels = pirna_genes_tb$other_name))

    # return
    return(fpkm_filt)
    
  }) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(animal = factor(animal, levels = c("mouse", "rat", "golden_hamster", "cow", "bull_testis")))

# plot
pirna_genes_plot <- 
  ggplot() +
  geom_bar(data = pirna_genes_fpkm, 
           mapping = aes(x = other_name, y = mean, fill = animal),
           width = 0.8, stat = "identity",
           position = position_dodge(width = 0.8)) +
  geom_errorbar(data = pirna_genes_fpkm,
                mapping = aes(x = other_name, ymin = mean - SEM, ymax = mean + SEM, color = animal), 
                width = 0.4, size = 1, 
                position = position_dodge(width = 0.8)) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  ggsave(filename = file.path(outpath, "piRNA_pathway_genes.FPKM.animals.barplot.png"), width = 15, height = 10)


### save FPKM values
# long table to wide
fpkm_wide <- 
  pirna_genes_fpkm %>% 
  dplyr::select(gene_id, gene_name, other_name, animal, mean) %>% 
  tidyr::pivot_wider(id_cols = c("gene_name", "other_name"), names_from = animal, values_from = mean) %>% 
  dplyr::arrange(gene_name)

# save
readr::write_csv(fpkm_wide, path = file.path(outpath, "piRNA_proteins.FPKM.csv"))
readr::write_csv(pirna_genes_fpkm, path = file.path(outpath, "piRNA_proteins.FPKM.long.csv"))

# SEM
pirna_genes_fpkm %>% 
  dplyr::select(gene_id, gene_name, other_name, animal, SEM) %>% 
  tidyr::pivot_wider(id_cols = c("gene_name", "other_name"), names_from = animal, values_from = SEM) %>% 
  dplyr::arrange(gene_name) %T>% 
  readr::write_csv(., path = file.path(outpath, "piRNA_proteins.SEM.csv"))