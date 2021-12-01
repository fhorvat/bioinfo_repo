### INFO: 
### DATE: Mon Apr 27 22:19:10 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
options(width = 185)

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/Analysis/2020_paper/small_RNA_reads_annotation/vfranke_test")

######################################################## LIBRARIES
library(data.table)
library(stringr)
library(ggplot2)
library(doMC)
library(GenomicAlignments)
library(reshape2)
library(edgeR)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
Assigner <- function(`.path`, `.name`){
  if(! is.character(`.path`) | !is.character(`.name`))
    stop('Both arguments should be characters!')
  
  load(`.path`)
  assign(`.name`, get(ls()[1]), parent.frame())
}

GffToGRanges <- function(gff, filter=NULL){
  
  library(plyr)
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

BamName <- function(bam.file){
  sub('.bam','',basename(bam.file))
}

AnnotateRanges <- function(r1, l, ignore.strand=FALSE,type = 'precedence', null.fact = 'None',collapse.char=':'){
  
  if(! class(r1) == 'GRanges')
    stop('Ranges to be annotated need to be GRanges')
  
  if(! all(sapply(l, class) == 'GRanges'))
    stop('Annotating ranges need to be GRanges')
  
  if(!type %in% c('precedence','all'))
    stop('type may only be precedence and all')
  
  require(data.table)
  require(GenomicRanges)
  cat('Overlapping...\n')
  if(class(l) != 'GRangesList')
    l = GRangesList(lapply(l, function(x){values(x)=NULL;x}))
  a = suppressWarnings(data.table(as.matrix(findOverlaps(r1, l, ignore.strand=ignore.strand))))
  a$id = names(l)[a$subjectHits]
  a$precedence = match(a$id,names(l))[a$subjectHits]
  a = a[order(a$precedence)]
  
  if(type == 'precedence'){
    cat('precedence...\n')
    a = a[!duplicated(a$queryHits)]
    annot = rep(null.fact, length(r1))
    annot[a$queryHits] = a$id
  }
  if(type == 'all'){
    cat('all...\n')
    a = a[,list(id=paste(unique(id),collapse=collapse.char)),by='queryHits']
    annot = rep(null.fact, length(r1))
    annot[a$queryHits] = a$id
    
  }
  return(annot)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "."

# gtf path
gene.path <- "/common/DB/vfranke/Base/GenomeAnnotation/mm9/Genes_ensembl/150523_mm9_Mus_musculus.NCBIM37.67.gtf"

# rmsk path
reps.path <- "/common/DB/vfranke/Base/GenomeAnnotation/mm9/RepeatMasker/mm9.RepeatMasker.RData"

# miRBase path
mirbase.path <- "/common/DB/vfranke/Base/MirBase/mm/mm9.LO.mm10.mirbase.form.gff"

# get bam file path
bamfile <- "./oocytesmallRNA_19to30.cleaned.fq.bam"

######################################################## READ DATA
# read repeatMasker
Assigner(reps.path, 'reps')

# read genes
genes <- data.table::fread(gene.path)
genes <-  genes[!str_detect(genes$V1, 'NT')]
genes <-  genes[!str_detect(genes$V1, 'MT')]
genes$V1 <- paste('chr', genes$V1, sep = '')
genes$V9 <- str_replace_all(genes$V9, '"', '')
gff <- GffToGRanges(as.data.frame(genes), 'exon')

# read miRBase
mirbase <- read.table(mirbase.path, header = FALSE, sep = '\t', stringsAsFactors = F)
mirbase <- mirbase[mirbase$V3 == "miRNA", ]
mirbase <- na.omit(mirbase)
mirbase <- GffToGRanges(mirbase)

######################################################## MAIN CODE
### get all the features from genomic .gtf
# exons
exon <- gff[gff$gene_biotype == 'protein_coding']

# introns
intron <- unlist(range(split(exon, exon$transcript_id)))

# pseudogenes
pseudo <- gff[gff$gene_biotype == 'pseudogene']

# exons on another strand
exon.anti <- exon
strand(exon.anti) <- ifelse(strand(exon.anti) == '+', '-', '+')

# lincRNA 
linc.exon <- gff[gff$gene_biotype == 'lincRNA']
linc.intron <- unlist(range(split(linc.exon, exon$transcript_id)))
linc.exon.anti <- linc.exon
strand(linc.exon.anti) <- ifelse(strand(linc.exon.anti) == '+', '-', '+')

# rRNA
rRNA.ens <- gff[gff$gene_biotype == 'rRNA']
rRNA.rep <- reps[reps$repClass == 'rRNA']

# tRNA
tRNA <- reps[reps$repClass == 'tRNA']

# creates list of annotations 
lreps <- list(MT = reps[str_detect(reps$repName, 'MT[ABCD]')],
              SINE = reps[reps$repClass == 'SINE'],
              LINE = reps[reps$repClass == 'LINE'],
              LTR = reps[reps$repClass == 'LTR'])
lannot1 <- list(rRNA.ens = rRNA.ens,
                rRNA.repmask = rRNA.rep,
                tRNA = tRNA,
                mirna = mirbase,
                pseudo = pseudo)
lannot2 <- list(exon.anti = exon.anti,
                linc.exon.anti = linc.exon.anti,
                exon = exon,
                linc.exon = linc.exon)
lannot <- c(lannot1, lreps, lannot2)

# create GRangesList (remove mcols values with lapply call)
lannot <- GRangesList(lapply(lannot, function(x){values(x) = NULL; x}))


### annotate reads
# get .bam name
name <- BamName(bamfile)

# read all reads, include number of mismatches tag
reads <- readGAlignments(bamfile, use.names = TRUE)

# remove reads mapping to mitochondria, get GRanges
seqlevels(reads, pruning.mode = "coarse") <- setdiff(seqlevels(reads), 'chrM')
g <- granges(reads, use.mcols = FALSE)
g$name <- paste(name, names(g), sep = '.')

# annotates ranges
gannot <- AnnotateRanges(g, lannot)

# creates data.table of annotations
tab <- data.table(name = names(reads))

# counts number of alignments for each read
tab <- tab[, .N, by = name]
tab$width <- width(reads)[match(tab$name, names(reads))]
tab$short <- tab$width >= 21 & tab$width <= 23
tab$uniq <- ifelse(tab$N == 1, 'Uniq', 'Mult')
tab <- tab[tab$N <= 5]

# add annotation, filter
tab$annot <- gannot[match(tab$name, names(reads))]
tab$annot <- factor(tab$annot, levels = c(names(lannot), 'None'), ordered = TRUE)
tabm <- tab[order(tab$annot), ]
tabm <- tabm[!duplicated(tabm$name)]

# summarize, cast to wide
dtabm <- tabm[tabm$short,length(unique(name)), by = list(annot, uniq)]
dtabm <- dcast.data.table(uniq~annot, data = dtabm)

# save
readr::write_delim(x = dtabm, path = file.path(outpath, str_c("read_class", name, "txt", sep = ".")), delim = "\t")
