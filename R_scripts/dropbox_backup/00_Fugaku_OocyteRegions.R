### INFO: R Script
### DATE: 05.08.2014
### AUTHOR: Vedran Franke
### PATH: /home/members/vfranke/Projects/PSvoboda_Oocyte_Transcriptome/Scripts/GenomeActivation/ExonAnalysis/CountReadsExonIntron.R
rm(list = ls()); gc()

################################################################################### WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/Fugaku/")

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
inpath <- "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2"

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/Fugaku/output/documentation"

################################################################################### SCRIPT PARAMS
# register workers for parallel computation
registerDoMC(21)

################################################################################### MAIN
files <- list.files(inpath, full.names = T, recursive = T, pattern = "bam$")
files <- files[(str_detect(files, "GV") | str_detect(files, "MII")) & !str_detect(files, "uniq")]

# create coverage in MII stages (.WE and .PA)
l.bam = list()
for(i in 2:length(files)){
  
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
  
  l.bam[[name]] <- regs
}

# add coverage together
l.regs <- lapply(l.bam, function(x) Reduce("+", x))

# take regions which have coverage of at least 5 and reduce them 
s <- lapply(l.regs, function(x) GenomicRanges::reduce(as(IRanges::slice(x, lower = 5), "GRanges")))

# reduce two stages to one 
d <- unlist(GRangesList(s))
d <- reduce(d)

# write as .bed 
write.table(as.data.frame(d)[, 1:3], 
            file.path(outpath, paste("OocyteRegions.lower.", i, ".bed", sep = "")), 
            row.names = F, col.names = F, quote = F, sep = "\t")








