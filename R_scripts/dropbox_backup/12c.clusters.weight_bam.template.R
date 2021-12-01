### INFO: filter bam, save as RData
### DATE: Fri Aug 31 14:24:21 2018
### AUTHOR: vfranke, modified by Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Analysis/smallRNA_clusters")

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
bam_name <- basename(bam_path) %>% str_remove_all(., ".bam")

# alignment count path
read_alignments_path <- file.path(inpath, str_c(bam_name, ".read_counts.txt.gz"))
  
# get chromosomes and their lengths from bam
chrs <- chrFinder(bam_path)

######################################################## READ DATA
# read alignment counts
read_alignments <- fread(read_alignments_path, col.names = c("N", "qname"))
  
######################################################## MAIN CODE
### add weighted NH values to reads
# register workers for parallel processing 
registerDoMC(10)

# read in parallel alignments from bam file, add values
l.bam <- foreach(chr = chrs$chr) %dopar% {
  
  # read bam, filter, get GRanges
  print(chr)
  which <- GRanges(chr, IRanges(1, chrs[chrs$chr == chr, ]$chrlen))
  param <- ScanBamParam(which = which, what = c("qname"))
  bam <- readGAlignments(bam_path, param = param)
  g <- granges(bam)
  
  # add weighted and raw nh value to multimapping reads
  if(length(g) > 0){
    values(g)$nh <- round((1 / read_alignments[match(values(bam)$qname, qname), N]), 3)
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

# save as .RData file
saveRDS(r, file = file.path(outpath, stringr::str_c(bam_name, ".weighted_reads.RDS")))

