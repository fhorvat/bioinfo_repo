### INFO: 
### DATE: Sat Sep 01 13:48:06 2018
### AUTHOR: vfranke (/home/members/vfranke/Projects/PSvoboda_Dicer/Scripts), modified by Filip Horvat
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

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

lib.path <- "/home/members/vfranke/MyLib/RFun"
source(file.path(lib.path, 'BamWorkers.R'))
source(file.path(lib.path, 'FileLoader.R'))
source(file.path(lib.path, 'FormatConverters.R'))
source(file.path(lib.path, 'ChipSeq.Functions.R'))
source(file.path(lib.path, 'ScanLib.R'))

######################################################## FUNCTIONS
### reads gff and converts it to GRanges (with optional filtering)
GffToGRanges <- function(gff, filter = NULL){
  
  if(ncol(gff) != 9)
    stop("Number of columns does not match gff format")
  
  if(any(gff[,5] < gff[,4])){
    warning("gff file contains ranges with negative widths...")
    gff = gff[gff[,5] > gff[,4],]
  }
  
  if(!is.null(filter)){
    if(filter %in% gff[,3]){
      cat("Filtering", filter, "features...\n")
      gff = gff[gff[,3] == filter,]
    }else{
      stop("The given feature is not present in the gff file")
    }
  }
  
  
  cat('Getting the feature ids...\n')
  s = strsplit(gff$V9, split=';')
  z = sapply(s, length)
  a = split(s, z)
  gff = gff[order(z),]
  l = lapply(a, function(x){
    d = sub('^ ','', unlist(x, use.names=F))
    d = sub('^.+? ','',d)
    m = matrix(d, ncol = length(x[[1]]), byrow=T)
    colnames(m) = sub(' .+$','',sub('^ ','', x[[1]]))
    m})
  ids = rbind.fill(lapply(l, data.frame))
  
  cat('Constructing the granges...\n')
  gff$V7[!gff$V7 %in% c('+','-')] = '*'
  granges = GRanges(seqnames = gff[,1],
                    IRanges(gff[,4],gff[,5]),
                    strand = gff[,7],
                    frame = gff[,8],
                    feature.type = gff[,3],
                    .id = 1:nrow(gff))
  
  values(granges) = cbind(values(granges), DataFrame(ids)[granges$.id,])
  values(granges)$.id = NULL
  return(granges)
}


### converts bed files to GenomicRanges
BedToGRanges <- function(bed, values = FALSE, seqlen = NULL){
  
  # order bed
  bed <- bed[order(bed[, 1], bed[, 2]), ]
  
  # check if there is strand information in bed
  if(sum(str_detect(names(bed), "strand")) == 1){
    
    strand <- as.character(bed[, str_detect(names(bed), "strand")])
    
    if(!all(strand %in% c("+", "-" , "*"))){
      stop('unallowed strand character found')
    }
    
  }else{
    
    strand <- "*"
    
  }
  
  # create GRanges
  ranges <- GRanges(seqnames = as.character(bed[, 1]),
                    ranges = IRanges(start = as.numeric(bed[, 2]), end = as.numeric(bed[, 3])),
                    strand = strand,
                    .id = 1:nrow(bed))
  
  # set seqlength
  if(!is.null(seqlen)){
    
    gseqlev <- GRanges(names(seqlen), IRanges(1, seqlen))
    ranges <- subsetByOverlaps(ranges, gseqlev, type = "within")
    seqlevels(ranges) <- names(seqlen)
    seqlengths(ranges) <- seqlen
    
  }
  
  # set values 
  if(values == TRUE && ncol(bed) > 3){
    
    col.ids <- setdiff(names(bed), c( "seqnames", "ranges", "strand", "seqlengths", "start", "end", "width",  "element", "chr"))
    values(ranges)[, col.ids] <- bed[values(ranges)$".id", col.ids]
    
  }
  
  values(ranges)$'.id' <- NULL
  
  # return
  return(ranges)
  
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
inpath <- "/common/WORK/vfranke/Projects/Dicer_RNASeq_29062012/Data/Mapped/ShrimpStrata"

# set outpath
outpath <- getwd()

# trimmed path
cut.path <- "/common/WORK/vfranke/Projects/Dicer_RNASeq_29062012/Data/Cutadapt_Trimmed2"

# rdata path
rdata <- file.path(outpath, "RData4")
dir.create(rdata, showWarnings = F)

# width path
wpath <- file.path(outpath, "Width")
dir.create(wpath, showWarnings = F)

# annotation path
apath <- file.path(outpath, "Annotation")
dir.create(apath, showWarnings = F)

# rdata stats path
stats.rdata.path <- file.path(outpath, "RDataStats")
dir.create(stats.rdata.path, showWarnings = F)

# ensembl genes path
ens.genes.path <- "/common/DB/vfranke/Base/GenomeAnnotation/mm9/Genes_ensembl/25012012.EnsemblGenesAll_BioMart.form.txt"

# exons path
exons.path <- "/common/DB/vfranke/Base/GenomeAnnotation/mm9/Genes_ensembl/22012012_EnsemblGenes_BioMart.form.gtf"

# miRNAs path
mirna.path <- "/common/DB/vfranke/Base/MirBase/mm/mm9.LO.mm10.mirbase.form.gff"

# repeatMasker path
repeat.path <- "/common/DB/vfranke/Base/GenomeAnnotation/mm9/RepeatMasker/mm9.RepeatMasker.RData"

# bam files path
bam.files <- list.files(inpath, pattern = 'bam$', full.names = T, recursive = T)
bam.files <- bam.files[str_detect(bam.files, "s_")]

######################################################## READ DATA
# microRNAs table
mirna <- read.table(mirna.path, header = F, stringsAsFactors = F)
mirna[, 9] <- str_replace_all(mirna[, 9], "=", " ")
mirna <- mirna[!(is.na(mirna[, 4]) | is.na(mirna[, 5])) ,]
mirna <- GffToGRanges(mirna)

# exons table
exons <- read.table(exons.path, header = F, sep = "\t", stringsAsFactors = F)
exons <- GffToGRanges(exons)

# genes table
genes <- read.table(ens.genes.path, header = T, sep = "\t")
genes <- BedToGRanges(genes, values = T)
genes <- genes[!values(genes)$biotype %in% c("pseudogene", "miRNA")]
exons <- exons[values(exons)$gene_id %in% values(genes)$ens.gene.id]
gr <- reduce(genes)
ex <- reduce(exons)

# repeats table
Assigner(repeat.path, "reps")
reps <- updateObject(reps)
reps <- reps[!str_detect(as.vector(seqnames(reps)), "random")]
seqlevels(reps) <- unique(as.character(seqnames(reps)))
reps <- reduce(reps)

######################################################## MAIN CODE
# set cuts for read lengts
cuts <- c(7, 15, 18, 20, 23, 34, 35, 49)

# add annotation to list
l.annot <- list(miRNA = mirna,
                Genes = ex,
                Repeats = reps)

# clean memory
rm(reps, gr, genes, exons, ex); gc()

# loop thorugh bams
for(i in 1:length(bam.files)){
  
  bam.file <- bam.files[i]
  name <- BamName(bam.file)
  print(name)
  
  # get chromosomes and their lengths from bam
  chrs <- chrFinder(bam.file)

  # classify reads from bam in parallel 
  registerDoMC(2)

  l <- foreach(chr = chrs$chr) %dopar% {
    
    # read alignments
    print(chr)
    which.ranges <- GRanges(chr, IRanges(1, chrs[chrs$chr == chr, ]$chrlen))
    param <- ScanBamParam(what = c("qname"), which = which.ranges, flag = scanBamFlag(isUnmappedQuery = FALSE), tag = c("NM"))
    bam <- readGAlignments(bam.file, param = param)
    
    # filter reads with more than 2 mismatches
    bam <- bam[values(bam)$NM <= 2]
    
    # annotate 
    d <- do.call(cbind, lapply(l.annot, function(x) suppressWarnings(countOverlaps(bam, x))))
    d[d > 1] <- 1
    d <- t(t(d) * 1:ncol(d))
    d[d == 0] <- max(d)+1
    da <- t(data.frame(d))
    mi <- do.call(pmin, split(da, 1:ncol(d)))
    ind <- c(names(l.annot), 'Other')[mi]
    dt <- data.table(width = qwidth(bam), qname = values(bam)$qname, ind = ind, NM = values(bam)$NM)
    rm(ws, s, bam, d, da, s, mi, ind); gc()
    
    # return
    return(dt)
    
  }
  
  
  # bind data.tables in list to one data.table
  r <- rbindlist(l)
  
  # assign read lengths to cuts 
  tab <- r[, .N, by = qname]
  m <- match(r$qname, tab$qname)
  r$nh <- tab$N[m]
  r$cut <- cut(r$width, cuts, include.lowest = T)
  
  # save as RData
  save(r, file = file.path(rdata, paste(name, "data.table.RData", sep = ".")))
  
}



