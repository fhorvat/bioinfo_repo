### INFO: classifies reads (not hierarchically) across repeatMasker
### DATE: Sun Sep 09 00:43:15 2018
### AUTHOR: vfranke, modified by Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/rmsk_counts/%EXPERIMENT")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(Rsamtools)
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

# genome path
genome_dir <- "/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/genomes/%EXPERIMENT"

# repeatMasker table path
exons_path <- list.files(path = genome_dir, pattern = "rmsk.*[0-9]{6}.clean.fa.out.gz", full.names = T)

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/Mapped/%EXPERIMENT"

# bam files path
bam_paths <- list.files(path = mapped_path, pattern = "*.total.bam$", full.names = T)


######################################################## READ DATA
# read repeatMasker
rmsk_df <- readr::read_delim(rmsk_path, delim = "\t")

######################################################## MAIN CODE
# prepare repeatMasker
rmsk_gr <- 
  rmsk_df %>% 
  dplyr::mutate(repFamily = ifelse(is.na(repFamily), repClass, repFamily)) %>% 
  GenomicRanges::GRanges(.) 

# set number of workers for parallel counting
registerDoMC(10)

### loop through all bam files
for(bam_path in bam_paths){
  
  # get bam name
  bam_name <- basename(bam_path) %>% str_remove_all(., ".bam")
  
  # get chromosomes and their lengths from bam
  chrs <- chrFinder(bam_path)
  
  # annotate reads in parallel 
  reads_class <- foreach(chr = chrs$chr) %dopar% {
    
    # read alignments
    print(chr)
    which.ranges <- GRanges(chr, IRanges(1, chrs[chrs$chr == chr, ]$chrlen))
    param <- ScanBamParam(what = c("qname"), which = which.ranges, flag = scanBamFlag(isUnmappedQuery = FALSE))
    bam <- readGAlignments(bam_path, param = param)
    
    # find overlaps with repeatMasker
    hits <- findOverlaps(bam, rmsk_gr, ignore.strand = T)
    
    # get names of reads and repClass/repFamily
    hits_df <- 
      tibble(qname = values(bam[queryHits(hits)])$qname, 
             repFamily = values(rmsk_gr[subjectHits(hits)])$repFamily, 
             repClass = values(rmsk_gr[subjectHits(hits)])$repClass, 
             repName = values(rmsk_gr[subjectHits(hits)])$repName) %>% 
      dplyr::filter(!duplicated(qname))
    
    # return
    return(hits_df)
    
  }
  ?findOverlaps
  # bind tibbles in list to one tibble
  reads_class_df <- 
    dplyr::bind_rows(reads_class) %>% 
    dplyr::filter(!duplicated(qname))
  
  # remove from memory
  rm(reads_class); gc()
  
  # summarize by repClass
  sum_repClass <-
    reads_class_df %>%
    dplyr::group_by(repClass) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(sample_id = bam_name)
  
  # summarize by repFamily
  sum_repFamily <-
    reads_class_df %>%
    dplyr::group_by(repFamily) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(sample_id = bam_name)
  
  # summarize by repName
  sum_repName <-
    reads_class_df %>%
    dplyr::group_by(repName) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(sample_id = bam_name)
  
  # save as RData
  saveRDS(object = list(repClass = sum_repClass, repFamily = sum_repFamily, repName = sum_repName), 
          file = file.path(outpath, str_c("rmsk.read_class", bam_name, "RDS", sep = ".")))
  
}



