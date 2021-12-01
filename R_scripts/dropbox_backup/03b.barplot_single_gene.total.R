### INFO: 
### DATE: Thu Apr 25 16:16:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/developmental_profile_expression/plots")

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

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 93


### genome
# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)


### experiment
# experiment list
experiment_list <- c("Fugaku", 
                     "Veselovska_2015_GenomeBiol_GSE70116")

# set base experiment path
base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/developmental_profile_expression")

# documentation path
documentation_path <- file.path(base_path, "Documentation")

# summarizedExperiment path
se_main_path <- file.path(base_path, "summarizedExperiments")

# FPKM path
fpkm_main_path <- file.path(base_path, "FPKMs")


### documentation and data
# sample table path
sample_table_path <- list.files(documentation_path, ".*sampleTable.csv", full.names = T)

# fpkm paths
fpkm_paths <- 
  list.files(path = fpkm_main_path, pattern = str_c("ensembl.", ensembl_version, ".*\\.FPKM_statistics\\.csv"), full.names = T) %>% 
  .[str_detect(., str_c(experiment_list, collapse = "|"))]


######################################################## READ DATA
# read sample table
sample_table <- data.table::fread(sample_table_path)

# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read counts
fpkm_tb <- 
  purrr::map(fpkm_paths, function(path){
    
    # read and transform to long format, join with sample table, calculate mean, SD and standard error
    data.table::fread(path) %>% 
      .[, experiment := str_extract(path, str_c(experiment_list, collapse = "|"))] %>% 
      .[]
    
  }) %>% 
  rbindlist(.) 

######################################################## MAIN CODE
# clean stage names
stage_clean <- 
  tibble(stage = c("nonGrowing_oocytes", "growing_oocytes_d8_14", "growing_oocytes_d15", "GV_oocytes",
                   "GV.WE.PE", "MII.WE.PE", "1cell.WE.PE", "2cell.WE.PE", "4cell.WE.PE", "Blast.WE.PE", "Molura.WE.PE"), 
         stage_clean = c("non-growing oocyte", "growing oocyte days 8-14", "growing oocyte day 15", "GV ", 
                         "GV", "MII", "1-cell", "2-cell", "4-cell", "Blastula", "Morula")) %>%
  mutate(stage_clean = factor(stage_clean, levels = stage_clean))

### filter by gene ID (lnc1)
# set gene id
chosen_gene_id <- "ENSMUSG00000078365"

# filter table
fpkm_tb_filt <- 
  fpkm_tb[gene_id == chosen_gene_id] %>% 
  as.tibble(.) %>% 
  right_join(., stage_clean, by = "stage") %>% 
  mutate_if(is.numeric, funs(ifelse(is.na(.), 0, .))) %>% 
  dplyr::left_join(., genes_info, by = "gene_id")

# write table
readr::write_csv(fpkm_tb_filt,
                 file.path(outpath, str_c("developmental_profile.barplot.FPKM", chosen_gene_id, "total_RNA.csv", sep = ".")))

### split to Veselovska and Fugaku
purrr::map(c("Fugaku", "Veselovska"), function(experiment_name){
  
  # visualize
  ggplot(fpkm_tb_filt %>% dplyr::filter(str_detect(experiment, experiment_name)), aes(x = stage_clean, y = avg_fpkm, fill = stage)) + 
    geom_bar(stat = "identity") +
    # geom_errorbar(aes(ymin = avg_fpkm - SD, ymax = avg_fpkm + SD), width = 0.2) +
    # geom_errorbar(aes(ymin = avg_fpkm - SE, ymax = avg_fpkm + SE), width = 0.2) +
    # scale_fill_grey() +
    scale_fill_viridis_d() +
    ylab("FPKM") + 
    xlab("") +
    # ggtitle(str_c(chosen_gene_id, " developmental profile")) +
    guides(fill = FALSE) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
          axis.text.y = element_text(size = 10, vjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggsave(filename = file.path(outpath, str_c("developmental_profile.barplot.FPKM", chosen_gene_id, "total_RNA", experiment_name, "viridis.png", sep = ".")), width = 8, height = 10)
  
})


