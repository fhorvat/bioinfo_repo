### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/CriGriPICR.LINE1_consensus.oocyte_RNA.STAR")

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

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq"

# documentation path
documentation_path <- file.path(base_path, "Data/Documentation")

# mapped path
mapped_path <- file.path(base_path, "Data/Mapped/STAR_mesAur1")

# sample table path
sample_tb_path <- list.files(documentation_path, pattern = "*.sampleTable.csv", full.names = T)

# stats and tracks path
stats_tb_path <- list.files(mapped_path, pattern = "*.stats_and_tracks.csv", full.names = T)

# deviate table path
deviate_tb_path <- 
  list.files(inpath, "*bam\\.L*", full.names = T) %>% 
  .[!str_detect(., "F66|F53|bai")]

######################################################## READ DATA
# read data
deviate_tb_full <- purrr::map(deviate_tb_path, function(path){
  
  # load full dataframe with fread and name the columns - 20 cols
  frame <- 
    fread(cmd = paste('grep -v ^#', path), data.table = FALSE, header = FALSE,
          colClasses = list(character = c(1, 2, 4, 14, 15, 18, 19, 20))) %>% 
    as_tibble(.) %>% 
    set_names(., c("TEfam", "sample_id", "pos", "refbase", "A", "C", "G", "T",
                   "coverage", "phys_coverage", "hq_coverage", "snp", "refsnp", "int_del", "int_del_freq",
                   "truncation_left", "truncation_right", "insertion", "deletion", "annotation")) %>% 
    dplyr::mutate(sample_id = sample_id %>% basename(.) %>% str_remove(., "\\.fastq\\.fused\\.sort\\.bam|\\.24to31nt\\.bam|\\.bam")) %>% 
    dplyr::select(-c("int_del", "int_del_freq", "truncation_left", "truncation_right", "insertion", "deletion", "annotation"))
  
  # return
  return(frame)
  
}) %>% 
  dplyr::bind_rows(.)

# read sample table
sample_tb <- readr::read_csv(sample_tb_path)

# read stats table
stats_tb <- readr::read_csv(stats_tb_path)

######################################################## MAIN CODE
# set color palette
base_colors <- c("A" = "#1b9e77", "C" = "#d95f02", 
                 "G" = "#7570b3", "T" = "#e6ab02",
                 "in" = "#66a61e", "del" = "#e7298a")

# create table with genotype and library size
sample_tb_filt <-
  sample_tb %>% 
  dplyr::select(sample_id, genotype) %>% 
  dplyr::left_join(., stats_tb %>% dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA), by = "sample_id") %>% 
  dplyr::mutate(library_size = (library_size / 10e5)) %>%
  dplyr::mutate(genotype = str_remove(genotype, "Mov10l_"),
                genotype = factor(genotype, levels = c("WT", "HET", "KO")))

# get unique repeat names
repeat_names <- unique(deviate_tb_full$TEfam)

### get data for plot 
purrr::map(repeat_names, function(repeat_name){
  
  # filter table
  deviate_tb <- 
    deviate_tb_full %>% 
    dplyr::filter(TEfam == repeat_name)
  
  # limit L1-1_CGr family plot
  if(repeat_name == "L1-1_CGr"){
    deviate_tb %<>%
      dplyr::filter(pos < 6000)
  }
  
  # get SNPs
  snps_tb <- 
    deviate_tb %>% 
    dplyr::filter(refsnp | snp) %>% 
    dplyr::select(sample_id, repName = TEfam, pos, refbase, A:T) %>% 
    dplyr::mutate(pos = pos + 1) %>% 
    tidyr::pivot_longer(A:T, names_to = "base_col", values_to = "counts_col") %>% 
    dplyr::filter(base_col != refbase) %>% 
    dplyr::select(sample_id, repName, base_col, pos_col = pos, counts_col)  %>% 
    dplyr::left_join(., sample_tb_filt, by = "sample_id") %>%
    dplyr::mutate(counts_col = round((counts_col / library_size), 3)) %>%
    dplyr::ungroup(.) %>% 
    dplyr::mutate(sample_id = str_remove_all(sample_id, "s_GV_Mov10l_|\\_r[0-9]\\.SE")) %>% 
    dplyr::mutate(sample_id = factor(sample_id, levels = c("WT_F70", "WT_So777", 
                                                           "HET_F58", "HET_F63", 
                                                           "KO_F52", "KO_F59", "KO_F77")))
  
  # add genotype and library size to table
  coverage_tb <-
    deviate_tb %>% 
    dplyr::left_join(., sample_tb_filt, by = "sample_id") %>% 
    dplyr::select(sample_id, repName = TEfam, pos, coverage, hq_coverage, genotype, library_size) %>% 
    dplyr::mutate(coverage = round((coverage / library_size), 3),
                  hq_coverage = round((hq_coverage / library_size), 3)) %>%
    dplyr::ungroup(.) %>% 
    dplyr::mutate(sample_id = str_remove_all(sample_id, "s_GV_Mov10l_|\\_r[0-9]\\.SE")) %>% 
    dplyr::mutate(sample_id = factor(sample_id, levels = c("WT_F70", "WT_So777", 
                                                           "HET_F58", "HET_F63", 
                                                           "KO_F52", "KO_F59", "KO_F77")))
  
  # creates base plot object
  famlen <- max(coverage_tb$pos) + 1
  lim <- ceiling(famlen / 50)
  
  # polygon total coverage
  pg <- tibble(px = coverage_tb$pos, 
               py = coverage_tb$coverage, 
               genotype = coverage_tb$genotype, 
               sample_id = coverage_tb$sample_id) 
  
  # polygon hq coverage
  pg_hq <- tibble(px = coverage_tb$pos, 
                  py = coverage_tb$hq_coverage, 
                  genotype = coverage_tb$genotype, 
                  sample_id = coverage_tb$sample_id) 
  
  # create plot
  baseplot <- 
    ggplot() +
    geom_rect(data = pg, aes(xmin = px, xmax = px + 1, ymin = 0, ymax = py), fill = "lightgrey") +
    geom_rect(data = pg_hq, aes(xmin = px, xmax = px + 1, ymin = 0, ymax = py), fill = "grey") +
    geom_bar(data = snps_tb, aes(x = pos_col, y = counts_col, fill = base_col), stat = "identity", width = 1.5) +
    facet_grid(rows = vars(sample_id)) +
    ylab("") +
    xlab("") +
    ggtitle(repeat_name) +
    scale_x_continuous(limits = c(-lim, famlen + lim)) +
    scale_fill_manual(values = base_colors, breaks = c("A", "C", "G", "T", "in", "del")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  # save plot
  ggsave(plot = baseplot, 
         filename = file.path(outpath, str_c("LINE1", repeat_name, "CriGriPICR.5K_to_7k.ORFs.consensus", 
                                             "Mov10l1", "oocyte_RNA", "individual", "png", sep = ".")), 
         width = 15, height = 7.5)
  
  # return 
  return(repeat_name)
  
})

