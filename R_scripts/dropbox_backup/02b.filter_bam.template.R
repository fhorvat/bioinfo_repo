### INFO: filter bam, save as RData
### DATE: Tue Jan 07 15:04:53 2020
### AUTHOR: vfranke, modified by Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd(".")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(doMC)
library(data.table)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
### assign RData file to variable
Assigner <- function(`.path`, `.name`){
  
  if(!is.character(`.path`) | !is.character(`.name`)){
    stop('Both arguments should be characters!')
  }
  
  load(`.path`)
  
  assign(`.name`, get(ls()[1]), parent.frame())
  
}


### scans chromosomes in bam file
chrFinder <- function(bam.path, filter = FALSE, output = "data.frame"){
  
  # scan bam
  s <- scanBamHeader(bam.path)
  st <- s[[1]]$text
  st <- do.call(rbind, st[names(st) == "@SQ"])
  st[, 1] <- str_replace(st[,1], "SN:", "")
  st[, 2] <- str_replace(st[,2], "LN:", "")
  
  # filter
  if(filter == TRUE){
    st <- st[!str_detect(st[, 1], "random")]
  }
  
  # output 
  if(output == 'data.frame'){
    
    vst <- data.frame(chr = st[, 1], chrlen = as.numeric(st[, 2]), stringsAsFactors = F)
    
  }else{
    
    vst <- as.numeric(st[, 2])
    names(vst) <- st[, 1]
    
  }
  
  # return
  return(vst)
  
}


######################################################## PATH VARIABLES
### in and out
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# arguments from command line
experiment <- args$experiment
threads <- as.numeric(args$threads)
genome_path <- args$genome_path
info_path <- args$info_path
exons_path <- args$exons_path
rmsk_path <- args$rmsk_path
mirbase_path <- args$mirbase_path
bam_path <- args$bam_path
bam_name <- args$bam_name
stats_path <- args$stats_path
nm_limit <- as.numeric(args$nm_limit)

######################################################## READ DATA
# load stats data from .RData files
Assigner(stats_path, "rdata")

######################################################## MAIN CODE
# set cut-off values
# nh.lim <- 4

# get chromosomes and their lengths from bam
chrs <- chrFinder(bam_path)

### save filtered bam files into a RData format
# subset reads - 21-23nt, number of mismatches, number of multimappers
rdata.sub <- subset(rdata, subset = (cut == "(20,23]" & NM <= nm_limit))
nm <- rdata.sub[, .N, by = qname]
m <- nm$N[match(rdata.sub$qname, nm$qname)]
rdata.sub$nh <- m
# rdata.sub <- subset(rdata.sub, subset = (nh <= nh.lim))

# get names of reads which passed filtering
qname.all <- rdata[, length(unique(qname[NM == nm_limit]))]
qname <- unique(rdata.sub$qname)

# clear memory
rm(rdata, nm, m); gc()

# read in parallel alignments from bam file which passed filtering
l.stats <- list()
registerDoMC(threads)

l.bam <- foreach(chr = chrs$chr) %dopar% {
  
  # read bam, filter, get GRanges
  print(chr)
  l.counts <- list()
  which <- GRanges(chr, IRanges(1, chrs[chrs$chr == chr, ]$chrlen))
  param <- ScanBamParam(which = which, what = c("qname"), tag = "NM")
  bam <- readGAlignments(bam_path, param = param)
  bam <- bam[values(bam)$qname %in% qname]
  g <- granges(bam)
  
  # add weighted and raw nh value to multimapping reads
  if(length(g) > 0){
    values(g)$nh <- (1 / rdata.sub[, nh[match(values(bam)$qname, qname)]])
    values(g)$nh_raw <- 1
  }
  
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
save(l, file = file.path(outpath, stringr::str_c(bam_name, "filt.RData", sep = ".")))

