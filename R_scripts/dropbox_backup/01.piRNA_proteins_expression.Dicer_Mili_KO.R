### INFO: 
### DATE: Sun Oct 06 17:21:37 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/datasets/2019_Sep/Analysis/expression/piRNA_pathway_genes")

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
inpath <- file.path(getwd(), "..")

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# FPKM table path
fpkm_path <- list.files(inpath, ".*\\.FPKM_long\\.csv", full.names = T)

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read FPKM table
fpkm_tb <- 
  readr::read_csv(fpkm_path) %>% 
  dplyr::left_join(., genes_info, by = "gene_id")

######################################################## MAIN CODE
# Aravin et al. (10.1371/journal.pgen.1000764): 
# MILI-TDRD1 - The first type of granules, pi-bodies, contains the MILI-TDRD1 module of the piRNA pathway
# MIWI2-TDRD9-MAEL - The second type of granules, piP-bodies, harbors the MIWI2-TDRD9-MAEL module of the piRNA 
# MVH - mouse VASA homolog (MVH) protein, an RNA helicase
# 
# Ernst et al. (10.1038/s41467-017-01049-7): 
# Maelstrom (MAEL) - piRNA precursor molecules are then exported into the cytoplasm, which might be facilitated by Maelstrom (MAEL)
# PLD6 (MITOPLD) - the endonuclease that generates the 5' end of primary piRNAs = Zuccchini
# Tdrkh - The process is facilitated by Tdrkh, a Tudor and KH domain-containing protein, whose absence causes the accumulation of 31–37 nucleotide (nt) long intermediates
# HENMT1 - Mature 3' ends are then 2'-O methylated by the RNA methyltransferase HENMT1
# RNF17 - The disengagement of MIWI from the ping-pong cycle was shown to be due to active repression of the ping-pong cycle by RNF17, a Tudor family member, in meiotic cells
# A-MYB - Transcription of these clusters is facilitated by the ancestral transcription factor A-MYB
# DDX4 - Vasa (also known as Ddx4) is an ATP-dependent RNA helicase

# # Pandey et al. (10.1073/pnas.1316316110)
# TDRD12 - Tudor domain containing 12 (TDRD12) is essential for secondary PIWI interacting RNA biogenesis in mice
# 
# Vagin et al. (10.1101/gad.1814809): 
# Piwi proteins were found in complex with PRMT5/WDR77, an enzyme that dimethylates arginine residues
# 
# Chen et al. (10.1073/pnas.0911640106): 
# Tdrd1, Tdrd4/RNF17, and Tdrd6, Tdrkh/Tdrd2
# In this analysis, we observed several previously known Piwi-interacting proteins or piRNA pathway components such as Ddx4, Kif17, Fxr1, and Mov10l1

# Siomi et al. (10.1101/gad.1899210)
# PIWI proteins have been shown recently to contain symmetrical dimethyl arginines (sDMAs), 
# and this modification is mediated by the methyltransferase PRMT5 (also known as Dart5 or Capsuleen)

### piRNA-pathway associated proteins
# create table
pirna_genes_tb <- 
  tibble(gene_name = c("Piwil2", "Piwil1", "Piwil4", 
                       "Mov10l1", "Ddx4", "Mael", 
                       "Pld6", "Henmt1",
                       "Tdrd1", "Tdrkh", "Rnf17", "Tdrd6", "Tdrd9", "Tdrd12", 
                       "Prmt5", "Wdr77", "Kif17"),
         other_name = c("MILI", "MIWI", "MIWI2", 
                        "MOV10L1", "DDX4/VASE", "MAEL/Maelstorm", 
                        "PLD6/MITOPOLD/Zucchini", "HENMT1", 
                        "TDRD1", "TDRD2/TDRKH", "TDRD4/RNF17", "TDRD6", "TDRD9", "TDRD12", 
                        "PRMT5/Dart5/Capsuleen", "WDR77", "KIF17"))


### FPKM table
# get FPKM of chosen genes
pirna_genes_fpkm <- 
  fpkm_tb %>% 
  dplyr::filter(gene_name %in% pirna_genes_tb$gene_name) %>% 
  dplyr::filter(sample_id != "s_GV_MILI_r1.PE") %>% 
  dplyr::mutate(genotype = str_remove_all(sample_id, "_old|_r[0-9]+.PE|s_GV_")) 

# data statistics
pirna_genes_fpkm_stats <- 
  pirna_genes_fpkm %>% 
  group_by(gene_name, genotype) %>% 
  summarize(mean = mean(fpkm), 
            sd = sd(fpkm), 
            length = length(fpkm)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(SEM = sd / sqrt(length)) %>% 
  dplyr::left_join(., pirna_genes_tb, by = "gene_name") %>% 
  dplyr::mutate(gene_name = factor(gene_name, levels = pirna_genes_tb$gene_name),
                other_name = factor(other_name, levels = pirna_genes_tb$other_name),
                genotype = factor(genotype, levels = c("WT", "MILI", "SOM", "DBL")))

# plot
pirna_genes_plot <- 
  ggplot() +
  geom_bar(data = pirna_genes_fpkm_stats, 
           mapping = aes(x = other_name, y = mean, fill = genotype),
           width = 0.8, stat = "identity",
           position = position_dodge(width = 0.8)) +
  geom_errorbar(data = pirna_genes_fpkm_stats,
                mapping = aes(x = other_name, ymin = mean - SEM, ymax = mean + SEM, color = genotype), 
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
  ggsave(filename = file.path(outpath, "piRNA_pathway_genes.FPKM.Dicer_Mili_KO.barplot.png"), width = 10, height = 10)

