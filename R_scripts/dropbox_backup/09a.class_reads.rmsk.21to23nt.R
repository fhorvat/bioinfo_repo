### INFO: classifies reads (not hierarchically) across repeatMasker
### DATE: Sun Sep 09 00:43:15 2018
### AUTHOR: vfranke, modified by Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/smallRNA_clusters")

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
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

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
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# get repeatMasker path
rmsk_path <- list.files(genome_path, pattern = "rmsk.*[0-9]{6}.clean.fa.out.gz", full.names = T)

# bam files path
# bam_path <- "%BAM_PATH"
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/T3T_DcrTrans_2011/Data/Mapped/Shrimp_mm10/s_T3T_DcrO_tran_r1.SE.mis_0.bam"
bam_name <- basename(bam_path) %>% str_remove_all(., ".bam")
experiment_name <- str_remove(bam_path, "/Data.*") %>% basename(.)

######################################################## READ DATA
# read repeatMasker
rmsk_df <- readr::read_delim(rmsk_path, delim = "\t")

######################################################## MAIN CODE
# prepare repeatMasker
rmsk_gr <- 
  rmsk_df %>% 
  dplyr::mutate(repFamily = ifelse(is.na(repFamily), repClass, repFamily)) %>% 
  GenomicRanges::GRanges(.) 

# get chromosomes and their lengths from bam
chrs <- chrFinder(bam_path)

# annotate reads in parallel 
registerDoMC(1)

reads_class <- foreach(chr = chrs$chr) %dopar% {
  
  chr <- "chr1"
  
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
           repClass = values(rmsk_gr[subjectHits(hits)])$repClass) %>% 
    dplyr::filter(!duplicated(qname))
  
  # return
  return(hits_df)
  
}

# bind tibbles in list to one tibble
reads_class_df <- 
  dplyr::bind_rows(reads_class) %>% 
  dplyr::filter(!duplicated(qname))

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

# save as RData
saveRDS(object = list(repClass = sum_repClass, repFamily = sum_repFamily), 
        file = file.path(outpath, str_c("rmsk.read_class", bam_name, "RDS", sep = ".")))



