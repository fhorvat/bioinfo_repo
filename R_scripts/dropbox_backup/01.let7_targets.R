### INFO: 
### DATE: Tue Sep 03 21:58:11 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/let7_targets")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# targetScan path
targetscan_dir <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/accessory_data/targetScan"

# get targetScan tables paths
targetscan_path <- file.path(targetscan_dir, c("targetScan.mouse.7.1.20181109.conserved_family.conserved_and_unconserved_targets.mouse.txt.gz", 
                                               "targetScan.mouse.7.1.20181109.conserved_family.conserved_targets.mouse.txt.gz"))

# Fugaku's FPKM path
fugaku_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Fugaku/Analysis/expression/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.Fugaku.FPKM.csv"
  
# CNOT6L FPKM path
cnot6l_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/expression/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.CNOT6L.avg.FPKM.csv"

######################################################## READ DATA
# # read genes info
# genes_info <- readr::read_csv(genes_info_path)

# read targetScan table
targetscan_list <- 
  purrr::map(targetscan_path, function(path) readr::read_delim(path, delim = "\t")) %>% 
  set_names(c("all", "conserved"))

# read Fugaku's FPKM values
fugaku_tb <- readr::read_csv(fugaku_path) 

# read CNOT6L FPKM values
cnot6l_tb <- readr::read_csv(cnot6l_path)

######################################################## MAIN CODE
### get let7 targets, save
# open workbook
wb <- openxlsx::createWorkbook()

# separately get and write all targets and only conserved targets
purrr::map(names(targetscan_list), function(conservation){
  
  # get targets, join with FPKM tables
  targetscan_tb <- 
    targetscan_list[[conservation]] %>% 
    dplyr::filter(str_detect(`miR Family`, "let-7-5p")) %>% 
    purrr::set_names(., str_replace_all(colnames(.), " ", "_")) %>% 
    dplyr::select(gene_id = Gene_ID, mir_family = miR_Family, seed_match = Seed_match, PCT) %>% 
    dplyr::mutate(gene_id = str_remove(gene_id, "\\.[0-9]+$")) %>% 
    dplyr::left_join(., fugaku_tb %>% dplyr::select(gene_id, Fugaku_GV_FPKM = s_GV.WE.PE, gene_name:gene_description), by = "gene_id") %>% 
    dplyr::left_join(., cnot6l_tb %>% dplyr::filter(stage == "s_GV_WT") %>% dplyr::select(gene_id, CNOT6L_GV_FPKM = avg_fpkm), by = "gene_id") %>% 
    dplyr::select(gene_id, gene_name, Fugaku_GV_FPKM, CNOT6L_GV_FPKM, mir_family:PCT, coordinates:gene_description)
  
    # write data
    openxlsx::addWorksheet(wb = wb, sheetName = str_c(conservation, "_targets"))
    openxlsx::writeData(wb = wb, sheet = str_c(conservation, "_targets"), x = targetscan_tb)
    
    # return table
    return(targetscan_tb)

})

# write table
openxlsx::saveWorkbook(wb = wb, 
                       file = file.path(outpath, "let7_5p.oocyte_FPKM.ensembl_93.20180903.xlsx"), 
                       overwrite = TRUE)




