### INFO: 
### DATE: Thu Apr 08 19:22:10 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_rat/datasets/oocytes.Dicer_MT_HET.2021_Apr/Analysis/expression")

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

library(biomaRt)

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
genome_dir <- "/common/DB/genome_reference/rat/rn6.Rnor_6.0.GCA_000001895.4"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# # tissue expression path
# expression_path_1 <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_rat/datasets/tissues.Dicer_MT_HET.2021_Apr/Analysis/expression"
# expression_path_1 <- file.path(expression_path_1, "ensembl.99.Rnor_6.0.20200415.UCSCseqnames.FPKM.csv")

# oocyte expression path
expression_path_2 <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_rat/datasets/oocytes.Dicer_MT_HET.2021_Apr/Analysis/expression"
expression_path_2 <- file.path(expression_path_2, "ensembl.99.Rnor_6.0.20200415.UCSCseqnames.FPKM.csv")

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read expression tables
# expression_tb_1 <- readr::read_csv(expression_path_1)
expression_tb_2 <- readr::read_csv(expression_path_2)

######################################################## MAIN CODE
# ### get pseudogene annotation
# # set animal and Ensembl version
# ensembl_name <- "mmusculus"
# ensembl_release <- 99
# 
# ## Ensembl versions
# ensembl_url <-
#   tibble(ens_version = c(99, 93, 92, 91, 89, 86),
#          date = c("Jan 2020", "Jul 2018", "Apr 2018", "Dec 2017", "May 2017", "Oct 2016"),
#          URL_archive = c("http://jan2020.archive.ensembl.org", 
#                          "http://jul2018.archive.ensembl.org", 
#                          "http://apr2018.archive.ensembl.org",
#                          "http://dec2017.archive.ensembl.org",
#                          "http://may2017.archive.ensembl.org",
#                          "http://oct2016.archive.ensembl.org")) %>%
#   dplyr::filter(ens_version == ensembl_release) %$%
#   URL_archive
# 
# # load Mart of mouse database from ensembl
# # sometimes function useMart isn't able to connect to server and returns error, this chunck repeats useMart until there is no error
# mart <- "error"
# count <- 0
# while(class(mart) == "character"){
#   
#   count <- count + 1
#   print(str_c(ensembl_name, "_gene_ensembl", " ", count))
#   
#   # load ENSEMBL mart
#   mart <- tryCatch(expr = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = stringr::str_c(ensembl_name, "_gene_ensembl"), host = ensembl_url),
#                    error = function(x) return("error"))
#   
#   
#   # stop if count get too big
#   if(count > 2){
#     stop("Something's not right")
#   }
#   
# }
# 
# # get attributes table
# attr_tb <- 
#   listAttributes(mart) %>% 
#   as_tibble(.) %>% 
#   dplyr::filter(page == "feature_page")
# 
# # get info about genes
# ensembl_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "description", "source_name"), mart = mart)

# join tables
expression_tb <- 
  expression_tb_2 %>% 
  dplyr::select(gene_id, starts_with("s_")) %>% 
  tidyr::pivot_longer(data = ., cols = -gene_id, names_to = "sample_id", values_to = "fpkm") %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "MT_HET|WT"), 
                tissue = str_remove_all(sample_id, "^s_|_r1.SE") %>% 
                  str_remove(., genotype) %>% 
                  str_remove(., "_") %>% 
                  str_replace(., "GVDicer_", "oocyte")) %>% 
  dplyr::left_join(., genes_info, by = "gene_id")

# logFC table
logfc_tb <- 
  expression_tb %>% 
  dplyr::select(gene_id, genotype, fpkm, gene_biotype) %>% 
  tidyr::pivot_wider(id_cols = c(gene_id, gene_biotype), names_from = "genotype", values_from = "fpkm", names_prefix = "FPKM.") %>% 
  dplyr::mutate(mean = ((FPKM.MT_HET + FPKM.WT) / 2),
                lfc = log2(FPKM.MT_HET + 0.1) - log2(FPKM.WT + 0.1)) %>% 
  dplyr::mutate(gene_biotype = replace(gene_biotype, gene_biotype != "pseudogene", "other"), 
                gene_biotype = factor(gene_biotype, levels = c("pseudogene", "other"))) %>% 
  dplyr::arrange(desc(gene_biotype))

# write table
# readr::write_csv(logfc_tb, file = file.path(outpath, "ensembl.99.Rnor_6.0.20200415.oocyte.log2FC.csv"))

# PLOT
ma_plot <- 
  ggplot(data = logfc_tb, aes(x = mean, y = lfc, color = gene_biotype)) + 
  geom_point(size = 5, shape = 20) +
  scale_x_log10(limits = c(0.001, 1e5), 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_color_manual(values = c("pseudogene" = "red3", "other" = "black")) +
  # scale_y_continuous(breaks = c(-5:3)) +
  # guides(color = FALSE) +
  xlab("average expression") +
  ylab("log2(Dicer MT HET / Dicer WT)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2), 
        axis.title.y = element_text(size = 15, vjust = 0.3), 
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, "MA_plot.ensembl.99.Rnor_6.0.20200415.oocyte.log2FC.png"),
       plot = ma_plot, width = 15, height = 10)

