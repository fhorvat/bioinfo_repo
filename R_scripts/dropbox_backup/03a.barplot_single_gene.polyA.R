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
experiment_list <- c("Deng_2014_Science_GSE45719", 
                     "Smallwood_2011_NatGenet_PRJEB2547",
                     "Yamaguchi_2013_CellRes_GSE41908", 
                     "Gan_2013_NatCommun_GSE35005", 
                     "ENCODE_2014_Nature_GSE49417")

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

# # read ENSEMBL reduced exons
# exons_gr <- readRDS(file = exons_path)

######################################################## MAIN CODE
# clean stage names
stage_clean <- 
  tibble(stage = c("PGC_9.5", "PGC_11.5", "PGC_13.5_m", "PGC_13.5_f", 
                   "primitive_SG_A", "SG_B", "leptotene_SC", "pachytene_SC", "round_ST", 
                   "d10_oocyte", "FG_GV_oocyte", 
                   "zygote", "2C", "4C", "8C", "16C", "mid_blast",
                   "placenta_adult8wks"), 
         stage_clean = c("PGC 9.5", "PGC 11.5", "male PGC 13.5", "female PGC 13.5", 
                         "primitive SG A", "SG B", "leptotene SC", "pachytene SC", "round ST",
                         "small oocyte", "fully-grown oocyte", 
                         "zygote", "2-cell", "4-cell", "8-cell", "16-cell", "blastocyst 94h", 
                         "placenta")) %>%
  mutate(stage_clean = factor(stage_clean, levels = stage_clean))


### filter by gene ID (lnc1)
# set gene id
chosen_gene_id <- "ENSMUSG00000110001"

# filter table
fpkm_tb_filt <- 
  fpkm_tb[gene_id == chosen_gene_id] %>% 
  as.tibble(.) %>% 
  right_join(., stage_clean, by = "stage") %>% 
  mutate_if(is.numeric, funs(ifelse(is.na(.), 0, .)))

# write table
readr::write_csv(fpkm_tb_filt,
                 file.path(outpath, str_c("developmental_profile.barplot.FPKM", chosen_gene_id, "polyA_RNA.csv", sep = ".")))

# visualize
ggplot(fpkm_tb_filt, aes(x = stage_clean, y = avg_fpkm, fill = stage)) + 
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = avg_fpkm - SD, ymax = avg_fpkm + SD), width = 0.2) +
  # geom_errorbar(aes(ymin = avg_fpkm - SE, ymax = avg_fpkm + SE), width = 0.2) +
  # scale_fill_grey() +
  scale_fill_viridis_d() +
  ylab("avg. FPKM") + 
  xlab("") + 
  guides(fill = FALSE) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5), 
        axis.text.y = element_text(size = 10, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  ggsave(filename = file.path(outpath, str_c("developmental_profile.barplot.FPKM.SD", chosen_gene_id, "polyA_RNA", "viridis.png", sep = ".")), width = 15, height = 10)


