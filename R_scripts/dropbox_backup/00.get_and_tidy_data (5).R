### INFO: 
### DATE: Thu Jul 18 19:52:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/arrays/Joshi_2007_BMCDevBiol_GSE5558")

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

library(GEOquery)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

######################################################## READ DATA
### get data from GEO
# chip annotation
annotation <- getGEO("GPL3771", destdir = outpath, AnnotGPL = T)

# series matrix
eset <- getGEO("GSE5558", destdir = outpath)

######################################################## MAIN CODE
### series matrix from GEO includes only log2FC values and not individual array intensities,  
### so we need to get expression value for individual arrays
# get list of array expression values from GEO 
array_list <- 
  purrr::map(sampleNames(eset), getGEO, destdir = outpath) %>% 
  set_names(sampleNames(eset))

# get intensity values, change flag with criteria from paper, save as .csv
purrr::map(names(array_list), function(array_name){
  
  # replace whitespaces in column names
  array_list[[array_name]] %>% 
    Table(.) %>% 
    as_tibble(.) %>% 
    set_colnames(., str_replace_all(colnames(.), " ", "_") %>% 
                   str_replace_all(., "__", "_")) %>% 
    dplyr::mutate(Flags = ifelse(test = ((CH1_MEAN / CH1_BKD_Mean) < 1.5 |
                                           (CH2_MEAN / CH2_BKD_Mean) < 1.5 |
                                           CH1_MEAN < 200 | CH2_MEAN < 200),
                                 yes = -100, Flags)) %T>%
    readr::write_csv(., file.path(outpath, str_c(array_name, ".csv")))
  
})

# extract data about samples from list
sample_table <- 
  purrr::map(array_list, function(eset) eset@header) %>% 
  bind_rows(.) %>% 
  dplyr::select(geo_accession, source_name_ch1, source_name_ch2) %>% 
  dplyr::mutate(stage = str_extract(source_name_ch1, "E17.5|E12.5|E14.5|newborn"), 
                tissue = "ovary",
                genotype_ch1 = str_extract(source_name_ch1, "FIGLA|Wild-type"), 
                genotype_ch2 = str_extract(source_name_ch2, "FIGLA|Wild-type")) %>% 
  dplyr::mutate_at(.vars = vars(starts_with("genotype")), 
                   .funs = ~(str_replace_all(., c("FIGLA" = "Figla_KO", "Wild-type" = "WT")))) %>% 
  dplyr::select(geo_accession, tissue, stage, genotype_ch1, genotype_ch2, source_name_ch1, source_name_ch2) %T>% 
  readr::write_csv(., file.path(outpath, "Joshi_2007_BMCDevBiol_GSE5558.sampleTable.csv"))

# get connection between probe ID's and gene symbols
annotation_tidy <- 
  annotation %>% 
  Table(.) %>% 
  as_tibble(.) %>% 
  dplyr::select(ID, gene_symbol = `Gene symbol`, GenBank = `GenBank Accession`, gene_name = `Gene title`, coordinates = `Chromosome annotation`) %>% 
  tidyr::separate(gene_symbol, into = c("gene_symbol", "gene_symbol2"), sep = "///") %>% 
  tidyr::separate(coordinates, into = c("coordinates", "coordinates2"), sep = "///") %>% 
  dplyr::mutate_at(.vars = vars(starts_with("coordinates")), 
                   .funs = list(~str_replace(., "Chromosome ", "chr") %>% 
                                  str_replace(",.+ \\(", " ") %>% 
                                  str_replace("\\.\\.", " ") %>% 
                                  str_remove_all("\\)|, complement"))) %T>% 
  readr::write_csv(., file.path(outpath, "GPL3771.annot.tidy.csv"))

### download and save GenBank sequences 
# get GenBank sequences
genBank_sequences <- ape::read.GenBank(access.nb = unique(annotation_tidy$GenBank))

# save as .fasta
write.dna(genBank_sequences, 
          file = file.path(outpath, "GPL3771.annot.GenBank.fasta"), 
          format = "fasta", nbcol = 1, colw = 80)


