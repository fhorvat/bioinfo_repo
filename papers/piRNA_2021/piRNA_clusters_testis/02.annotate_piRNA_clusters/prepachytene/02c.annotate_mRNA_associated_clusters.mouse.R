### INFO: 
### DATE: Mon Aug 31 14:57:45 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/annotate_clusters")

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

library(openxlsx)
library(rtracklayer)
library(GenomicRanges)
library(biomaRt)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# clusters path
clusters_path <- file.path(inpath, "piRNA_associated_mRNAs.MesAur1.1k_windows.rpkm_cutoff.1_supplementary_pre-pachytene.20200831.FH.xlsx")

# genome path - hamster
genome_path <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"
gtf_path <- file.path(genome_path, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.gtf")
gene_info_path.hamster <- str_replace(gtf_path, "\\.gtf$", ".geneInfo.csv")

# genome path - mouse
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"
gtf_path <- file.path(genome_path, "ensembl.99.GRCm38.p6.20200415.UCSCseqnames.gtf.gz")
gene_info_path.mouse <- str_replace(gtf_path, "\\.gtf\\.gz$", ".geneInfo.csv")

# FPM tables path - spermatogonia from Garcia-Lopez, pachytenes from Papd7_Lacz_WT
fpm_GL_sense_path <- file.path(inpath, "ensembl.99.GRCm38.p6.20200415.UCSCseqnames.sense.GL.mouse.FPM.csv")
fpm_GL_antisense_path <- file.path(inpath, "ensembl.99.GRCm38.p6.20200415.UCSCseqnames.antisense.GL.mouse.FPM.csv")
fpm_Papd7_sense_path <- file.path(inpath, "ensembl.99.GRCm38.p6.20200415.UCSCseqnames.sense.Papd7.mouse.FPM.csv")
fpm_Papd7_antisense_path <- file.path(inpath, "ensembl.99.GRCm38.p6.20200415.UCSCseqnames.antisense.Papd7.mouse.FPM.csv")

######################################################## READ DATA
# read clusters table
clusters_tb <- 
  openxlsx::read.xlsx(clusters_path) %>% 
  as_tibble(.) %>% 
  unique(.)

# read geenes info
gene_info.hamster <- readr::read_csv(gene_info_path.hamster)
gene_info.mouse <- readr::read_csv(gene_info_path.mouse)

# read FPM tables
fpm_GL_sense <- readr::read_csv(fpm_GL_sense_path)
fpm_GL_antisense <- readr::read_csv(fpm_GL_antisense_path)
fpm_Papd7_sense <- readr::read_csv(fpm_Papd7_sense_path)
fpm_Papd7_antisense <- readr::read_csv(fpm_Papd7_antisense_path)

######################################################## MAIN CODE
### find hamster genes in mouse
# add by gene name
hamster_mouse_genes <- 
  clusters_tb %>% 
  dplyr::select(gene_id, gene_name) %>% 
  dplyr::left_join(., gene_info.mouse %>% dplyr::select(gene_name, gene_id.mouse = gene_id), by = "gene_name") %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::slice(1) %>% 
  dplyr::ungroup(.)

# set ensembl version
ensembl_version <- 99

# Ensembl versions
ensembl_url <-
  tibble(ens_version = c(100, 99, 98, 96, 95, 94, 93, 92, 91, 89, 86),
         date = c("Apr 2020", "Jan 2020", "Sep 2019", "Apr 2019", "Jan 2019", "Oct 2018", "Jul 2018", "Apr 2018", "Dec 2017", "May 2017", "Oct 2016"),
         URL_archive = c("www.ensembl.org",
                         "http://jan2020.archive.ensembl.org",
                         "http://sep2019.archive.ensembl.org",
                         "http://apr2019.archive.ensembl.org",
                         "http://jan2019.archive.ensembl.org",
                         "http://oct2018.archive.ensembl.org",
                         "http://jul2018.archive.ensembl.org",
                         "http://apr2018.archive.ensembl.org",
                         "http://dec2017.archive.ensembl.org",
                         "http://may2017.archive.ensembl.org",
                         "http://oct2016.archive.ensembl.org")) %>%
  dplyr::filter(ens_version == ensembl_version) %$%
  URL_archive

# load ENSEMBL mart
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mauratus_gene_ensembl", host = ensembl_url)

# get gene homologs
ensembl_homologs <-
  getBM(attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene"),
        filters = "ensembl_gene_id",
        values = hamster_mouse_genes$gene_id,
        mart = mart) %>%
  as_tibble(.) %>% 
  dplyr::filter(mmusculus_homolog_ensembl_gene != "") %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::slice(1) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::select(gene_id = ensembl_gene_id, gene_id.mouse_biomart = mmusculus_homolog_ensembl_gene)

# add to the table
hamster_mouse_genes %<>% 
  dplyr::left_join(., ensembl_homologs, by = "gene_id") %>% 
  dplyr::mutate(gene_id.mouse = ifelse(is.na(gene_id.mouse), gene_id.mouse_biomart, gene_id.mouse)) %>% 
  dplyr::select(gene_id, gene_id.mouse)


### add mouse expression to clusters table
# select samples
mouse_fpm <- 
  fpm_GL_sense %>% 
  dplyr::select(gene_id, mouse.spermatogonia_GL.sense = s_spermatogonia) %>% 
  dplyr::left_join(., fpm_GL_antisense %>% dplyr::select(gene_id, mouse.spermatogonia_GL.antisense = s_spermatogonia)) %>% 
  dplyr::left_join(., fpm_Papd7_sense %>% dplyr::select(gene_id, mouse.testis_Papd7_LacZ_WT.sense = s_testis_Papd7_LacZ_WT)) %>% 
  dplyr::left_join(., fpm_Papd7_antisense %>% dplyr::select(gene_id, mouse.testis_Papd7_LacZ_WT.antisense = s_testis_Papd7_LacZ_WT)) %>% 
  dplyr::rename(gene_id.mouse = gene_id)
  
# add mouse gene_id to the table
clusters_tb_mouse <- 
  clusters_tb %>% 
  dplyr::left_join(., hamster_mouse_genes, by = "gene_id") %>% 
  dplyr::left_join(., mouse_fpm, by = "gene_id.mouse")

# write
openxlsx::write.xlsx(clusters_tb_mouse, 
                     file = file.path(outpath, "piRNA_associated_mRNAs.MesAur1.1k_windows.rpkm_cutoff.1_supplementary_pre-pachytene.mouse_expression.20200903.FH.xlsx"))



