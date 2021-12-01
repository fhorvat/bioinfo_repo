### INFO: 
### DATE: Tue Jan 07 15:36:14 2020
### AUTHOR: vfranke, modified by Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2019_Dec/Analysis/expression/miRBase")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(data.table)
library(GenomicRanges)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)

######################################################## SOURCE FILES

######################################################## FUNCTIONS



######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gtf path
gtf_path <- file.path(genome_path, "miRBase.22.mm10.20181605.gff3")

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2019_Dec"

# FPM path
fpm_path <- file.path(inpath, "miRBase.22.mm10.20181605.FPM.csv")

# sample table path
sample_table_path <- file.path(base_path, "Data/Documentation")
sample_table_path <- list.files(sample_table_path, ".*\\.sampleTable\\.csv", full.names = T)

# library size path
library_size_path <- file.path(base_path, "Data/Mapped/STAR_mm10/7_perfect_reads") 
library_size_path <- list.files(library_size_path, "library_sizes.txt", full.names = T)

# pepa's annoatation path
pepa_path <- file.path(inpath, "mirAnnot.dt.rda")

######################################################## READ DATA
# read gtf info
mirbase_gr <- rtracklayer::import(gtf_path)

# read small RNA-seq expression
fpm_tb <- readr::read_csv(fpm_path)

# read sample table
sample_tb <- readr::read_csv(sample_table_path)

# load Pepa's annotation
load(pepa_path) 
  
######################################################## MAIN CODE
# tidy sample table
sample_tb %<>% 
  dplyr::select(sample_id, genotype)

# get miRBase annoation table
mirbase_tb <- 
  mirbase_gr %>% 
  as_tibble(.) %>% 
  tidyr::unite(coordinates, seqnames, start, end, sep = " ") %>% 
  dplyr::select(gene_id = Name, ID = ID, parent = Derives_from, type, coordinates) %>% 
  dplyr::filter(type == "miRNA") %>% 
  dplyr::left_join(., mirAnnot.dt %>% dplyr::select(ID, strand_type = type), by = "ID")
  
# get mean FPM per genotype
fpm_mean <- 
  fpm_tb %>% 
  tidyr::pivot_longer(-c(gene_id, coordinates), names_to = "sample_id", values_to = "fpm") %>% 
  dplyr::left_join(., sample_tb, by = "sample_id") %>% 
  dplyr::group_by(coordinates, genotype) %>% 
  dplyr::summarise(fpm = mean(fpm)) %>% 
  dplyr::ungroup(.)

# table for plot
plot_df <- 
  fpm_mean %>% 
  dplyr::filter(genotype != "HET") %>% 
  tidyr::pivot_wider(coordinates, names_from = "genotype", values_from = "fpm") %>% 
  dplyr::left_join(., mirbase_tb, by = "coordinates") %>% 
  dplyr::mutate(mean = ((WT + KO) / 2),
                lfc = log2(WT + 1) - log2(KO + 1))

# create ggplot object
ma_plot <- 
  ggplot() +
  geom_point(data = plot_df, aes(x = mean, y = lfc, color = strand_type), size = 3, shape = 20) +
  scale_x_log10(limits = c(0.01, 1e6), 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(limits = c(-6, 6)) +
  scale_color_manual(values = c("miRNA*" = "red", "miRNA" = "blue", "not_defined" = "grey50")) +
  # guides(color = T) +
  xlab("average expression (RPM)") +
  ylab("log2(DicerX WT / DicerX KO) E15.5") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save
ggsave(filename = file.path(outpath, str_c("miRBase", "MA_plot", "FPM", "DicerX_embryos", "png", sep = ".")), 
       plot = ma_plot,
       width = 10, height = 10)



