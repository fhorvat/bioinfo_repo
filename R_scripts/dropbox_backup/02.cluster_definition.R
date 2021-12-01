### INFO: 
### DATE: Fri Aug 31 15:10:23 2018
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

library(Cairo)
library(data.table)
library(GenomicRanges)
library(plyr)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS
### set date for outpath
OutpathDate <- function(name){
  paste(strsplit(as.character(Sys.time()), '\\s')[[1]][1], name, sep = '_')
}


### assign RData file to variable
Assigner <- function(`.path`, `.name`){
  
  if(!is.character(`.path`) | !is.character(`.name`)){
    stop('Both arguments should be characters!')
  }
  
  load(`.path`)
  
  assign(`.name`, get(ls()[1]), parent.frame())
  
}


### for each region sum weighted counts of all reads
GetOverlaps <- function(reg1, reg2, colname = NULL){
  
  if(is.null(colname)){
    stop("Colname needs to be defined")
  }
  
  fo <- data.table(as.matrix(findOverlaps(reg1, reg2)), ignore.strand=T)
  fo$weight <- values(reg2)[[colname]][fo$subjectHits]
  fo <- fo[, sum(weight), by = queryHits]
  v <- rep(0, length(reg1))
  v[fo$queryHits] <- fo$V1
  
  return(v)
  
}


### finds the clusters in short RNA data
FindClusters <- function(reads, total, rpkm1 = 5, rpkm2 = 10, clust.width = 100, e = 1e6){
  
  # reduce regions
  cat("Reducing the regions\n")
  dregs <- reduce(reads, ignore.strand = T)
  cat("clustnum:", length(dregs), "\n")
  
  # get weighted counts in each region
  co <- GetOverlaps(dregs, reads, "nh")
  
  # filtering based on RPKM
  cat("First step filtering\n")
  rpkm <- (co) * (e / total)
  rind <- rpkm > rpkm1
  wregs <- dregs[rind]
  values(wregs)$counts <- co[rind]
  values(wregs)$rpkm <- rpkm[rind]
  
  # join cluster closer than clust.width
  cat("Joining clusters\n")
  a.regs <- reduce(resize(wregs, width = width(wregs) + clust.width, fix = "center"), ignore.strand = T)
  a.regs <- resize(a.regs, width = width(a.regs) - clust.width, fix = "center")
  values(a.regs)$counts <- GetOverlaps(a.regs, reads, "nh")
  values(a.regs)$rpkm.c <- (values(a.regs)$counts) * (e / total)
  regs.sel <- a.regs
  
  # return 
  cat("Returning the clusters...\n")
  return(reg.sel = regs.sel)
  
}


### merges sets of cluster
ClusterMerge <- function(regs1, regs2, method = "intersect", ignore.strand = T){
  
  cat("method:", method, "\n")
  
  # intersect
  if(method == "intersect"){
    
    fo <- as.matrix(findOverlaps(regs1, regs2))
    o <- reduce(c(regs1[fo[, 1]], regs2[fo[, 2]]), ignore.strand = ignore.strand)
    
  }
  
  # union
  if(method == "union"){
    
    o <- reduce(c(regs1, regs2), ignore.strand = ignore.strand)
    
  }
  
  # return
  return(regs.sel = o)
  
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


### get annotation - .gtf to GenomicRanges
GetAnnotation <- function(l.paths){
  
  # get paths
  mirna.path <- l.paths$mirna.path
  exons.path <- l.paths$exons.path
  genes.path <- l.paths$ens.genes.path
  reps.path <- l.paths$repeat.path
  
  # read and tidy miRNA 
  mirna <- read.table(mirna.path, header = F, stringsAsFactors = F)
  mirna[, 9] <- str_replace_all(mirna[, 9], "=", " ")
  mirna <- mirna[!(is.na(mirna[, 4]) | is.na(mirna[, 5])), ]
  mirna <- GffToGRanges(mirna)
  
  # read exons
  exons <- read.table(exons.path, header = F, sep = "\t", stringsAsFactors = F)
  exons <- GffToGRanges(exons)
  
  # read genes
  genes <- read.table(genes.path, header = T, sep = "\t", stringsAsFactors = F)
  genes <- BedToGRanges(genes, values = T)
  genes <- genes[!values(genes)$biotype %in% c("pseudogene", "miRNA")]
  
  # filter exons
  exons <- exons[values(exons)$gene_id %in% values(genes)$ens.gene.id]
  
  # reduce genes and exons
  gr <- reduce(genes)
  ex <- reduce(exons)
  
  # reads repeatMasker
  Assigner(reps.path, "reps")
  reps <- updateObject(reps)
  reps <- reps[!str_detect(as.vector(seqnames(reps)), "random")]
  seqlevels(reps) <- unique(as.character(seqnames(reps)))
  reps <- reduce(reps)
  
  # return
  list(miRNA = mirna,
       Genes = ex,
       Repeats = reps)
  
}


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


### hierarchically annotate regions with list of annotations  
AnnotateRegions <- function(regs, l.annot){
  
  # count overlaps between regions and annotation
  d <- do.call(cbind, lapply(l.annot, function(x) suppressWarnings(countOverlaps(regs, x))))
  
  # set all regions which overlap more than 1 annotation of same type to 1
  d[d > 1] <- 1
  
  # multiply numbers in each column with index of that column
  d <- t(t(d) * 1:ncol(d))
  
  # set all regions which don't overlap with any annotation to highest number
  d[d == 0] <- max(d) + 1
  
  
  # for each region get lowest numbered annotation
  da <- t(data.frame(d))
  mi <- do.call(pmin, split(da, 1:ncol(d)))
  
  # get names of annotation for each region
  ind <- c(names(l.annot), "Other")[mi]
  
  # return
  return(ind)
  
}


### plot clusters in scatterplot
PlotClusters <- function(expr1, expr2, ind, outpath, name, xlab = "", ylab = "", mcw = NULL, cols = NULL){
  
  # set plotting parameters
  lwd <- 5
  cex.axis <- 2
  cex <- 1.5
  ind <- as.factor(ind)
  
  # set colors
  if(is.null(cols)){
    cols = c("red", "cornflowerblue", "black", "darkgray")
  }
  
  # plot 1
  outname.sq <- paste(name, "mcw", mcw, "sq.png", sep = ".")
  CairoPNG(file.path(outpath, outname.sq), width = 800, height = 800)
  plot(expr1, expr2, pch = 20, xlab = xlab, ylab = ylab, col = cols[as.numeric(ind)], axes = F, xlim = c(0, 6), ylim = c(0, 6), cex = cex)
  axis(1, at = 0:5, lwd = lwd, cex = cex.axis, tcl = -1, labels = F,  padj = 5)
  axis(2, at = 0:5, lwd = lwd, cex = cex.axis, tcl = -1, labels = F,  padj = 5)
  legend("topright", legend = levels(ind), fill = cols, bty = "n")
  dev.off()
  
}


######################################################## PATH VARIABLES
# set inpath 
inpath <- "/common/WORK/vfranke/Projects/Dicer_RNASeq_29062012/Results/ClusterDefinitionShrimp/Regions"

# set outpath
outdir <- getwd()

# create subdirs in outpath
outpath <- file.path(outdir, OutpathDate("Clust"))
dir.create(outpath, showWarnings = F)
plot.path <- file.path(outpath, "Plots")
dir.create(plot.path, showWarnings = F)
clust.path <- file.path(outpath, "Clust")
dir.create(clust.path, showWarnings = F)

# annotation paths
l.paths <- list(mirna.path = "/common/DB/vfranke/Base/MirBase/mm/mm9.LO.mm10.mirbase.form.gff",
                ens.genes.path = "/common/DB/vfranke/Base/GenomeAnnotation/mm9/Genes_ensembl/25012012.EnsemblGenesAll_BioMart.form.txt",
                repeat.path = "/common/DB/vfranke/Base/GenomeAnnotation/mm9/RepeatMasker/mm9.RepeatMasker.RData",
                exons.path = "/common/DB/vfranke/Base/GenomeAnnotation/mm9/Genes_ensembl/22012012_EnsemblGenes_BioMart.form.gtf")

# list count matrices RData files
l.files <- 
  list.files(inpath, full.names = T, pattern = "regions") %>% 
  .[str_detect(., "s_")] %>% 
  .[str_detect(., "nh")]

######################################################## READ DATA
# read count matrices
l.mat <- list()
for(i in 1:length(l.files)){
  
  l.file <- l.files[i]
  name <- str_replace(basename(l.file), ".nm.+$", "")
  print(name)
  Assigner(l.file, "l")
  l.mat[[name]] <- l
  
}

# read bams from count matrices
bams <- lapply(l.mat, function(x) x$reads)

# read annotation
l.annot <- GetAnnotation(l.paths)

######################################################## MAIN CODE
# set date
dat <- 
  str_replace_all(date(), "\\s+", "-") %>% 
  str_replace_all(., ":", "_")


### find clusters in bam files
# DicerO
rpkm1 <- 3
rpkm2 <- 3
width <- 50
o2.clust <- FindClusters(bams[["s_ES_MII.Dicer_Rsc_utra_r1"]], l.mat[["s_ES_MII.Dicer_Rsc_utra_r1"]]$total, rpkm1, rpkm2, width, 1e6)
o3.clust <- FindClusters(bams[["s_ES_MII.Dicer_Rsc_utra_r2"]], l.mat[["s_ES_MII.Dicer_Rsc_utra_r2"]]$total, rpkm1, rpkm2, width, 1e6)

# merge replicates
o.clust <- ClusterMerge(o2.clust, o3.clust)


### count reads in clusters
# for each bam count reads in clusters
d.cnts <- do.call(cbind, lapply(bams, function(x) GetOverlaps(o.clust, x, "nh")))

# get total values for each bam
total <- unlist(lapply(l.mat, "[[", "total"))

# normalize counts
d.cnts.norm <- data.frame(t(t(d.cnts) / total * 1e6))


### output DicerO table
# for each loci get percentage of reads mapping to +/- strand
cp <- GetOverlaps(o.clust, bams[["s_ES_MII.Dicer_Rsc_utra_r2"]][strand(bams[["s_ES_MII.Dicer_Rsc_utra_r2"]]) == "+"], "nh")
cm <- GetOverlaps(o.clust, bams[["s_ES_MII.Dicer_Rsc_utra_r2"]][strand(bams[["s_ES_MII.Dicer_Rsc_utra_r2"]]) == "-"], "nh")
perc <- cp / (cp + cm)

# create data.frame, set strand 
do <- data.frame(as.data.frame(o.clust), d.cnts.norm)
do$strand <- "*"
do$strand[perc > 0.8] <- "+"
do$strand[perc < 0.2] <- "-"

# annotate
do$class <- AnnotateRegions(o.clust, l.annot)

# save table
write.table(do, file.path(clust.path, str_c("O.combined", dat, "txt", sep = ".")), row.names = F, col.names = T, quote = F, sep = "\t")


### plot clusters
PlotClusters(expr1 = log10(do$s_ES_MII.Dicer_Rsc_utra_r1 + 1),
             expr2 = log10(do$s_ES_WT_tran + 1),
             ind = do$class,
             outpath = plot.path,
             name = "test",
             xlab = "DicerO-3",
             ylab = "DicerS-11")




