### INFO: classifies reads to RData data.table
### DATE:Tue Jan 07 14:11:08 2020
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

######################################################## READ DATA
# read info about genes
genes_info <- readr::read_csv(info_path)

# gtf exons by genes
exons_gr <- readRDS(file = exons_path)

# read repeatMasker
rmsk_df <- readr::read_delim(rmsk_path, delim = "\t")

# read miRBase gff
mirna_gr <- rtracklayer::import.gff(con = mirbase_path) 

######################################################## MAIN CODE
### prepare annotations - repeatMasker, exons of genes from ENSEMBL annotation, miRBase
# prepare exons
exons_gr %<>% 
  as.data.frame(.) %>% 
  tibble::as_tibble(.) %>% 
  dplyr::left_join(., genes_info %>% dplyr::select(gene_id, gene_biotype), by = "gene_id") %>% 
  dplyr::select(seqnames:strand, gene_id, gene_biotype) %>%
  dplyr::filter(!str_detect(gene_biotype, "pseudogene|miRNA")) %>% 
  GenomicRanges::GRanges(.)

# prepare repeats
rmsk_gr <- 
  rmsk_df %>% 
  GenomicRanges::GRanges(.) %>% 
  GenomicRanges::reduce(.)

# prepare mature miRNA
mirna_gr <- mirna_gr[mcols(mirna_gr)$type == "miRNA"]
mcols(mirna_gr) <- mcols(mirna_gr)[, c("Name")]
names(mcols(mirna_gr)) <- "gene_id"
mcols(mirna_gr)$gene_biotype <- "miRNA.mature"

# add annotation to list
l.annot <- list(miRNA = mirna_gr,
                Genes = exons_gr,
                Repeats = rmsk_gr)

# clean memory
rm(mirna_gr, exons_gr, rmsk_gr, rmsk_df); gc()


### annotate reads 
# set cuts for read lengts
cuts <- c(7, 15, 18, 20, 23, 34, 35, 49)

# get chromosomes and their lengths from bam
chrs <- chrFinder(bam_path)

# annotate reads in parallel 
registerDoMC(threads)

l <- foreach(chr = chrs$chr) %dopar% {
  
  # read alignments
  print(chr)
  which.ranges <- GRanges(chr, IRanges(1, chrs[chrs$chr == chr, ]$chrlen))
  param <- ScanBamParam(what = c("qname"), which = which.ranges, flag = scanBamFlag(isUnmappedQuery = FALSE), tag = c("nM"))
  bam <- readGAlignments(bam_path, param = param)
  
  # annotate 
  d <- do.call(cbind, lapply(l.annot, function(x) suppressWarnings(countOverlaps(bam, x))))
  d[d > 1] <- 1
  d <- t(t(d) * 1:ncol(d))
  d[d == 0] <- max(d)+1
  da <- t(data.frame(d))
  mi <- do.call(pmin, split(da, 1:ncol(d)))
  ind <- c(names(l.annot), "Other")[mi]
  dt <- data.table(width = qwidth(bam), qname = values(bam)$qname, ind = ind, NM = values(bam)$nM)
  rm(ws, s, bam, d, da, s, mi, ind); gc()
  
  # return
  return(dt)
  
}

# bind data.tables in list to one data.table
r <- data.table::rbindlist(l)

# count multimappers
tab <- r[, .N, by = qname]
m <- match(r$qname, tab$qname)
r$nh <- tab$N[m]

# assign read lengths to bins 
r$cut <- cut(r$width, cuts, include.lowest = T)

# save as RData
save(r, file = file.path(outpath, str_c(bam_name, "dt.RData", sep = ".")))
