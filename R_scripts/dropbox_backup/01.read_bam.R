### INFO: 
### DATE: Fri Aug 31 14:24:21 2018
### AUTHOR: vfranke, modified by Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/smallRNA_clusters/filtered_bams")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(doMC)
library(Cairo)
library(data.table)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

lib_path <- "/home/members/vfranke/MyLib/RFun"
source(file.path(lib_path, 'BamWorkers.R'))
source(file.path(lib_path, 'FormatConverters.R'))
source(file.path(lib_path, 'FileLoader.R'))
source(file.path(lib_path, 'FileLoader.R'))
source(file.path(lib_path, 'ScanLib.R'))

######################################################## FUNCTIONS
### gets name of bam file
BamName <- function(bam.file){
  sub(".bam", "", basename(bam.file))
}


### assign RData file to variable
Assigner <- function(`.path`, `.name`){
  
  if(!is.character(`.path`) | !is.character(`.name`)){
    stop('Both arguments should be characters!')
  }
  
  load(`.path`)
  
  assign(`.name`, get(ls()[1]), parent.frame())
  
}


######################################################## PATH VARIABLES
# set inpath 
inpath <- "/common/WORK/vfranke/Projects/Dicer_RNASeq_29062012/Data/Mapped/ShrimpStrata"

# set outpath
outpath <- getwd()

# bam path
bam.files <- list.files(inpath, pattern = "bam$", full.names = T, recursive = T)
bam.files <- bam.files[str_detect(bam.files, "s_")]

# list stats files
stats.path <- "/common/WORK/vfranke/Projects/Dicer_RNASeq_29062012/Results/MappingStats/ShrimpStrata/RData4"
stats.files <- list.files(stats.path, full.names = T, pattern = "RData")

######################################################## READ DATA

######################################################## MAIN CODE
# set cut-off values
nm.lim <- 0
nh.lim <- 4

# get lengths of chromosomes in genome
seqlen <- seqlengths(genome)
seqlen <- seqlen[!str_detect(names(seqlen), "_")]

### saves the filtered bam files into a RData format
l <- list()

for(i in 1:length(bam.files)){
  
  # prepare vectors
  l.samps <- list()
  bam.file <- bam.files[i]
  name <- BamName(bam.file)
  print(name)
  
  # load stats data from .RData files
  cat("Loading the rdata\n")
  rdata.file <- stats.files[str_detect(stats.files, name)]
  Assigner(rdata.file, "rdata")
  
  # subset reads - 21-23nt, number of mismatches, number of multimappers
  cat("Subsetting the rdata\n")
  rdata.sub <- subset(rdata, subset = (cut == "(20,23]" & NM <= nm.lim))
  nm <- rdata.sub[, .N, by = qname]
  m <- nm$N[match(rdata.sub$qname, nm$qname)]
  rdata.sub$nh <- m
  rdata.sub <- subset(rdata.sub, subset = (nh <= 4))
  
  # get names of reads which passed filtering
  cat("Getting the reads rdata\n")
  qname.all <- rdata[, length(unique(qname[NM == 0]))]
  qname <- unique(rdata.sub$qname)
  
  # clear memory
  rm(rdata, nm, m); gc()
  
  # read in parallel alignments from bam file which passed filtering
  l.stats = list()
  registerDoMC(2)
  
  l.bam <- foreach(chr = names(seqlen)) %dopar% {
    
    # read bam, filter, get GRanges
    print(chr)
    l.counts <- list()
    which <- GRanges(chr, IRanges(1, seqlen[chr]))
    param <- ScanBamParam(which = which, what = c("qname"), tag = "NM")
    bam <- readGAlignments(bam.file, param = param)
    bam <- bam[values(bam)$qname %in% qname]
    g <- granges(bam)
    
    # add weighted nh value to multimapping reads 
    values(g)$nh <- round(1 / rdata.sub[, nh[match(values(bam)$qname, qname)]], 2)
    
    # sanity check
    if(length(g) != length(bam)){
      stop('smth gone wrong!!!')
    }
    
    # return
    return(g)
    
  }
  
  # join all chromosomes 
  r <- unlist(GRangesList(l.bam))
  
  # create list with reads, total number of reads and number of 21-23nt reads
  l <- list(reads = r, total = qname.all, total.21 = length(qname))
  
  # save as .RData file
  save(l, file = file.path(outpath, paste(name, "nm", nm.lim, "nh", nh.lim, "regions.RData", sep = ".")))
  
  # clear memory
  rm(r, l.bam, l); gc()
  
}


