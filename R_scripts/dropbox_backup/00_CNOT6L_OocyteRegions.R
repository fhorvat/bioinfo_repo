### INFO: R Script
### DATE: 05.08.2014
### AUTHOR: Vedran Franke
### PATH: /home/members/vfranke/Projects/PSvoboda_Oocyte_Transcriptome/Scripts/GenomeActivation/ExonAnalysis/CountReadsExonIntron.R
rm(list = ls()); gc()

################################################################################### WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/")

################################################################################### LIBRARIES
lib.path <- "/common/WORK/fhorvat/R_library/vfranke"
source(file.path(lib.path, "FileLoader.R"))
source(file.path(lib.path, "FormatConverters.R"))
source(file.path(lib.path, "BamWorkers.R"))
source(file.path(lib.path, "ScanLib.R"))
library(data.table)
library(stringr)
library(doMC)
library(GenomicAlignments)
library(readr)
library(dplyr)
library(magrittr)

################################################################################### PATH VARIABLES
inpaths <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10"

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/documentation"

################################################################################### SCRIPT PARAMS
# register workers for parallel computation
registerDoMC(21)

################################################################################### TABLES
# CNOT6L experiment table
sample_table <- 
  read_csv("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/CNOT6L_sample_list_11919R_2015_10_29.csv", col_names = T) %>%
  dplyr::select(ID, stage = `Time Course`, treatment = `Treatment/Control`) %>%
  dplyr::mutate(name = str_c(ID, stage, treatment, sep = "_")) 

# library size
sample_table_path <- 
  tibble(sample_path = list.files(path = inpaths, pattern = "*.bam$", full.names = T, recursive = T)) %>%
  mutate(ID = str_replace_all(sample_path, "^/.*/|_.*", "")) %>%
  dplyr::select(ID, sample_path)

# join together, select only WT
sample_table <- 
  dplyr::left_join(sample_table, sample_table_path, by = "ID") %>% 
  dplyr::filter(treatment == "WT") %>% 
  dplyr::select(ID, name, sample_path)

################################################################################### MAIN
files <- 
  sample_table %>% 
  dplyr::filter(str_detect(name, "MII")) %$% 
  sample_path

# create coverage in MII stage
l_bam = list()
for(i in 1:length(files)){
  
  file = files[i]
  name = BamName(file)
  print(name)
  chrs = chrFinder(file)
  regs = foreach(chr = chrs$chr)%dopar%{
    
    print(chr)
    which.ranges <- GRanges(chr, IRanges(1, chrs$chrlen[chrs$chr == chr]))
    bam <- readGAlignments(file, param = ScanBamParam(which = which.ranges), use.names = TRUE)
    gbam <- granges(bam)
    return(coverage(gbam))
    
  }
  
  l_bam[[name]] <- regs
  
}

# add coverage together
l_regs <- lapply(l_bam, function(x) Reduce("+", x))

# take regions which have coverage of at least 5 and reduce them 
s <- lapply(l_regs, function(x) GenomicRanges::reduce(as(IRanges::slice(x, lower = 5), "GRanges")))

# reduce two stages to one 
d <- unlist(GRangesList(s))
d <- reduce(d)

# write as .bed 
write.table(as.data.frame(d)[, 1:3], 
            file.path(outpath, "CNOT6L_MII_OocyteRegions_lower5.bed"), 
            row.names = F, col.names = F, quote = F, sep = "\t")








