### INFO: filter bam, save as RData
### DATE: Fri Aug 31 14:24:21 2018
### AUTHOR: vfranke, modified by Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/smallRNA_clusters/T3T_DcrTrans_2011/filtered_bams")

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
# lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
# source(file.path(lib_path, "wideScreen.R"))
# source(file.path(lib_path, "headt.R"))
# source(file.path(lib_path, "asdf.R"))
# wideScreen()


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
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# bam path
bam_path <- "%BAM_PATH"
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/T3T_DcrTrans_2011/Data/Mapped/Shrimp_mm10/s_T3T_DcrO_tran_r1.SE.mis_0.bam"
bam_name <- basename(bam_path) %>% str_remove_all(., ".bam")

# list stats files
stats.path <- file.path(inpath, "../mapping_stats")
stats.files <- list.files(stats.path, full.names = T, pattern = "RData")

######################################################## READ DATA

######################################################## MAIN CODE
# set cut-off values
nm.lim <- 0
# nh.lim <- 4

# get chromosomes and their lengths from bam
chrs <- chrFinder(bam_path)

### save filtered bam files into a RData format
# load stats data from .RData files
rdata.file <- stats.files[str_detect(stats.files, bam_name)]
Assigner(rdata.file, "rdata")

# subset reads - 21-23nt, number of mismatches, number of multimappers
rdata.sub <- subset(rdata, subset = (cut == "(20,23]" & NM <= nm.lim))
nm <- rdata.sub[, .N, by = qname]
m <- nm$N[match(rdata.sub$qname, nm$qname)]
rdata.sub$nh <- m
# rdata.sub <- subset(rdata.sub, subset = (nh <= nh.lim))

# get names of reads which passed filtering
qname.all <- rdata[, length(unique(qname[NM == 0]))]
qname <- unique(rdata.sub$qname)

# clear memory
rm(rdata, nm, m); gc()

# read in parallel alignments from bam file which passed filtering
l.stats <- list()
registerDoMC(10)

l.bam <- foreach(chr = chrs$chr) %dopar% {
  
  # read bam, filter, get GRanges
  print(chr)
  l.counts <- list()
  which <- GRanges(chr, IRanges(1, chrs[chrs$chr == chr, ]$chrlen))
  param <- ScanBamParam(which = which, what = c("qname"), tag = "NM")
  bam <- readGAlignments(bam_path, param = param)
  bam <- bam[values(bam)$qname %in% qname]
  g <- granges(bam)
  values(g)$qname <- values(bam)$qname
  
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

  

