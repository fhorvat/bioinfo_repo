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
  fpkm_tb %>% 
  as_tibble(.) %>% 
  mutate_if(is.numeric, ~(ifelse(is.na(.), 0, .))) %>% 
  dplyr::filter(str_detect(experiment, "Veselovska")) %>% 
  dplyr::filter(stage %in% c("growing_oocytes_d8_14", "GV_oocytes")) %>% 
  dplyr::select(gene_id, stage, avg_fpkm) %>% 
  tidyr::pivot_wider(id_cols = gene_id, names_from = stage, values_from = avg_fpkm) %>% 
  dplyr::mutate(gene = ifelse(gene_id == chosen_gene_id, "Mos", "other")) %>% 
  dplyr::mutate(mean = ((GV_oocytes + growing_oocytes_d8_14) / 2),
                lfc = log2(GV_oocytes + 0.1) - log2(growing_oocytes_d8_14 + 0.1)) %>% 
  dplyr::left_join(., genes_info %>% dplyr::select(gene_id, gene_name, gene_biotype, gene_description), by = "gene_id") %>% 
  dplyr::filter(gene_biotype == "protein_coding")
# dplyr::select(gene_id, mean, lfc, gene)

# write table
fpkm_tb_filt %>%
  dplyr::select(-gene) %>% 
  readr::write_csv(.,
                   file.path(outpath, str_c("Veselovska.MA_plot.FPKM", "total_RNA.csv", sep = ".")))

# get Mos lfc
mos_lfc <- 
  fpkm_tb_filt %>% 
  dplyr::filter(gene_id == chosen_gene_id) %$% 
  lfc

# # normalize to Mos lfc
# fpkm_tb_filt %<>% 
#   dplyr::mutate(lfc = lfc - mos_lfc)

# result limits
results_limits <-
  fpkm_tb_filt %>% 
  dplyr::summarise(x_limit = mean %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.),
                   y_limit = lfc %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.))

# annotation table
annotations <- tibble(xpos = Inf,
                      ypos = mos_lfc,
                      annotateText = as.character(round(mos_lfc, 3)))

# plot
ma_plot <-
  ggplot() +
  geom_point(data = fpkm_tb_filt, aes(x = mean, y = lfc, color = gene, size = gene, alpha = gene), shape = 20) +
  geom_hline(yintercept = mos_lfc) +
  geom_text(data = annotations, aes(x = xpos, y = ypos, label = annotateText), 
            colour = "black", fontface = "italic", size = 2.5, 
            hjust = 1.2, vjust = -0.5) +
  scale_x_log10(limits = c(0.001, results_limits$x_limit),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(limits = c(-results_limits$y_limit, results_limits$y_limit),
                     breaks = c(-results_limits$y_limit:results_limits$y_limit)) +
  scale_colour_manual(values = c(Mos = "red2", other = "grey60")) +
  scale_alpha_manual(values = c(Mos = 1, other = 0.5)) +
  scale_size_manual(values = c(Mos = 5, other = 2.5)) +
  guides(color = guide_legend(override.aes = list(shape = 23, size = 5, fill = c("red2", "grey60"))),
         alpha = F,
         size = F) +
  xlab("mean expression") +
  ylab(str_c("log2 fold change: ", "GV oocyte", " / ", "growing oocyte (d8-d14)", "\n")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)) +
  theme(legend.title = element_blank())

# # turns off axis titles and legend
# ma_plot <-
#   ma_plot +
#   theme(legend.position = "none") +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank())

# save plot
ggsave(filename = file.path(outpath,
                            str_c("Veselovska.MA_plot.FPKM", chosen_gene_id, "total_RNA", "png", sep = ".")),
       plot = ma_plot, width = 12, height = 10)
