### INFO: 
### DATE: Sat Sep 01 13:48:06 2018
### AUTHOR: vfranke, modified by Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/smallRNA_clusters/mapping_stats")

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

# get gene info path
info_path <- list.files(genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.geneInfo.csv", full.names = T)

# get paths of reduced exons
exons_path <- list.files(genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.reducedExons.RDS", full.names = T)

# get repeatMasker path
rmsk_path <- list.files(genome_path, pattern = "rmsk.*[0-9]{6}.clean.fa.out.gz", full.names = T)

# get miRBase gff path
mirbase_path <- list.files(genome_path, pattern = "miRBase.*gff3", full.names = T)

# bam files path
bam_path <- "%BAM_PATH"
# bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/T3T_DcrTrans_2011/Data/Mapped/Shrimp_mm10/s_T3T_DcrO_tran_r1.SE.mis_0.bam"
bam_name <- basename(bam_path) %>% str_remove_all(., ".bam")
experiment_name <- str_remove(bam_path, "/Data.*") %>% basename(.)

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
  tibble::as.tibble(.) %>% 
  dplyr::left_join(., genes_info %>% dplyr::select(gene_id, gene_biotype), by = "gene_id") %>% 
  dplyr::select(seqnames:strand, gene_id, gene_biotype) %>%
  dplyr::filter(!str_detect(gene_biotype, "pseudogene|miRNA")) %>% 
  GenomicRanges::GRanges(.)

# prepare repeats
rmsk_gr <- 
  rmsk_df %>% 
  GenomicRanges::GRanges(.) %>% 
  reduce(.)

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
registerDoMC(2)

l <- foreach(chr = chrs$chr) %dopar% {
  
  # read alignments
  print(chr)
  which.ranges <- GRanges(chr, IRanges(1, chrs[chrs$chr == chr, ]$chrlen))
  param <- ScanBamParam(what = c("qname"), which = which.ranges, flag = scanBamFlag(isUnmappedQuery = FALSE), tag = c("NM"))
  bam <- readGAlignments(bam_path, param = param)
  
  # annotate 
  d <- do.call(cbind, lapply(l.annot, function(x) suppressWarnings(countOverlaps(bam, x))))
  d[d > 1] <- 1
  d <- t(t(d) * 1:ncol(d))
  d[d == 0] <- max(d)+1
  da <- t(data.frame(d))
  mi <- do.call(pmin, split(da, 1:ncol(d)))
  ind <- c(names(l.annot), "Other")[mi]
  dt <- data.table(width = qwidth(bam), qname = values(bam)$qname, ind = ind, NM = values(bam)$NM)
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
save(r, file = file.path(outpath, str_c(bam_name, "dt.RData2", sep = ".")))


