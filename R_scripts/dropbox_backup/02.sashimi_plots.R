### INFO: 
### DATE: Fri Jun 22 12:56:05 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/review/spliced_reads_sashimi")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicAlignments)
library(GenomicRanges)
library(grid)
library(rtracklayer)
library(gtable)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS
#' Use Bezier curves to make 'swoop' shapes represnting splice junctions
make.swoop <- function(df, scale = c(1.5e3, 15), ...) {
  
  # add ID value for each splice junction
  df$row <- seq_len(nrow(df))
  
  # calculate Bézier curves for each splice junction
  ddply(df, .(row), function(d) {
    
    # get width of splice junction
    w <- with(d, end - start)[1]
    
    # get middle point
    mid <- c(0.2, 0.8) * w + d$start[1]
    
    # calculate Bézier curves
    x <- c(d$start[1], mid, d$end[1])
    y <- c(0, -1 * rep(pmax(w / scale[1], scale[2]), 2), 0)
    rez <- as.data.frame(Hmisc::bezier(x, y))
    
    for (c in colnames(d)) {
      rez[, c] <- d[1, c]
    }
    
    return(rez)
    
  })
  
}

#' Render a GRanges of exons into a transcript ideogram
#' 
#' @param exons a \code{GRanges} containing exons; metadata column \code{"Parent"} groups them into distinct transcripts
#' @param at vertical offset for first transcript
#' @param height vertical size of each transcript
#' @param arrows if arrows desired along intron connectors (to show strand orientation), this should be a call
#' 	to \code{grid::arrow()} returning the desired sort of arrow
#' @param do.introns if \code{TRUE}, draw lines through introns to connect exons
#' @param colour line colour for intron connectors
#' @param colour.by name of metadata column to use for exon boxes
#' @param fill fill colour for exon boxes, if \code{colour.by} not specified
#' @param stroke line colour for borders of exon boxes
#' 
#' @value a list of \code{ggplot2} geoms, which can be added to an existing plot with \code{`+`}.
#' 
#' @details Designed for use with the output of \code{rtracklayer::import.gff3()} with GFFs from Ensembl.  For
#' 	other use cases, mileage may vary.
#' 

geom_tx <- function(exons, at = 0, height = 1, arrows = grid::arrow(length = unit(4, "pt"), type = "closed"),
                    do.introns = TRUE, colour.by = NULL, stroke = NA, colour = "black", fill = "black",
                    label = FALSE, label.size = 3, drop = TRUE, ...) {
  
  .make.tx.df <- function(ex) {
    
    if (!length(ex))
      return(NULL)
    
    if (drop)
      values(ex)$at <- as.numeric(ex$parent)*height + at
    else
      values(ex)$at <- as.numeric(ex$parent)
    values(ex)$height <- height
    if (drop)
      exx <- droplevels(as.data.frame(ex))
    else
      exx <- as.data.frame(ex)
    exx$panel <- "transcripts"
    if (is.null(colour.by) || is.na(colour.by)) {
      rez <- list( geom_rect(data = exx, aes(xmin = start, xmax = end, ymin = at-height*0.4, ymax = at+height*0.4),
                             fill = fill, colour = stroke, ...) )
    }
    else {
      exx$fill <- exx[ ,colour.by ]
      rez <- list( geom_rect(data = exx, aes(xmin = start, xmax = end, ymin = at-height*0.4, ymax = at+height*0.4, fill = fill),
                             colour = stroke, ...) )
    }
    
    if (label) {
      labels <- data.frame(xpos = max(exx$end) + 1000, ypos = exx$at[1]+exx$height[1]*0.4, label = exx$parent[1], panel = exx$panel[1])
      rez[[3]] <- geom_text(data = labels, aes(x = xpos, y = ypos, label = label), size = label.size, hjust = 0)
    }
    
    if (length(ex) > 1) {
      ## get exons from introns
      if (min(start(ex)) > 0)
        .introns <- GenomicRanges::gaps(ex)[-1]
      else
        .introns <-  GenomicRanges::gaps(ex)
      if (drop)
        introns <- droplevels(as.data.frame(.introns))
      else
        introns <- as.data.frame(.introns)
      #message(nrow(introns), " introns...")
      
      if (any(strand(ex) == "-")) {
        x <- introns$start
        introns$start <- introns$end
        introns$end <- x
      }
      
      #print(introns)
      introns$at <- ex$at[1]
      introns$height <- height
      introns$panel <- "transcripts"
      #if (drop)
      #	introns$middle <- with(introns, at+height*0.8/2)
      #else
      #	introns$middle <- with(introns, at + as.numeric(parent))
      rez <- c(rez, geom_segment(data = introns, aes(x = start, xend = end, y = at, yend = at),
                                 arrow = arrows, colour = colour, ...))
      
    }
    
    return(rez)
    
  }
  
  if (!inherits(exons, "GRanges"))
    exons <- makeGRangesFromDataFrame(exons, keep.extra.columns = TRUE)
  
  ## if input was GFF3, do everything on a per-tx basis
  if (!is.null(exons$Parent)) {
    
    if (is.null(exons$parent)) {
      ## rtracklayer::import.gff3() uses list format for Parent field; assume 1 parent, listed first
      if (inherits(exons$Parent,"List"))
        exons$parent <- sapply(exons$Parent, "[", 1)
      else
        exons$parent <- as.character(exons$parent)
    }
    else {
      ## respect existing 'parent' column
    }
    ## 'parent' needs to be a factor; force this to be so
    if (!is.factor(exons$parent))
      exons$parent <- reorder(factor(exons$parent), start(exons), min)
    
  }
  else {
    exons$parent <- 1
  }
  
  exons <- exons[ order(exons$parent) ]
  exl <- split(exons, exons$parent)
  lapply(exl, .make.tx.df)
  
}

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# get paths of genome reference files
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# get path of original bam file
bam_paths <- list.files(inpath, ".*CNOT6L.PE.total.bam")

# gtf path 
gtf_path <- file.path(genome_dir, "ensembl.86.GRCm38.p5.20180512.UCSCseqnames.gtf.gz")

######################################################## READ DATA
# read gtf
gtf <- import.gff(gtf_path)

######################################################## MAIN CODE
#' Make a 'sashimi plot': coverage plus splice events
#' 
#' @param aln a named list of \code{GAlignments} objects
#' @param tx a \code{GRanges} object, preferably as returned by \code{import.gff3()}
#' @param meta a dataframe of extra metadata; linked to samples by column \code{"iid"}, and column
#' 	\code{"reads"} is used to normalize coverage, if available
#' @param smooth integer; if >0, size of windows in which to calculate smoothed coverage estimate
#' @param min.coverage only show regions with coverage strictly greater than this
#' @param min.splices only show splice events with multiplicity strictly greater than this
#' @param max.coverage truncate coverage at this value
#' @param log.coverage logical; if \code{TRUE}, show coverage in log10 scale
#' @param splice.scale numeric vector of length 2 used for drawing splice-junction arcs; just play with it to find good values
#' @param colours a named vector of colours used for coverage plots; grey/black is default
#' 
#' @value a \code{grid} object with the completed plot

# CNOT6L coordinates
CNOT6L_gr <- GRanges(seqnames = "chr5", ranges = IRanges(start = 96070333, end = 96161547), strand = "-")

# deletion coordinates
del_gr <- tibble(seqnames = "chr5", start = 96075067, end = 96106332, strand = "-")

# read and subset bams
bam_GAlignment_list <- 
  purrr::map(bam_paths, function(bam_path){
    
    readGAlignmentsList(file = bam_path, param = ScanBamParam(what = "qname")) %>% 
      unlist(.) %>% 
      subsetByOverlaps(., ranges = CNOT6L_gr)
      
  }) %>% 
  set_names(., basename(bam_paths) %>% stringr::str_remove(., ".CNOT6L.PE.total.bam"))

aln <- bam_GAlignment_list
tx <- subsetByOverlaps(gtf, CNOT6L_gr)
meta <- data.frame(iid = basename(bam_paths) %>% stringr::str_remove(., ".CNOT6L.PE.total.bam"),
                   reads = c(18760601, 15465835, 15638192, 14707301, 18711816, 14716904))
min.coverage = 0
min.splices = 0
max.coverage = Inf
colours = NULL
colour.by = NULL
splice.scale = c(30000, 0.1)

aln <- bam_GAlignment_list[str_detect(names(bam_GAlignment_list), "KO")]
meta <-
  meta %>%
  dplyr::filter(str_detect(iid, "KO"))

## set boundaries of region to plot: approx. the range of all transcripts
zoom <- reduce(tx, min.gapwidth = 100e4)

## calculate raw coverage
cvg <- lapply(lapply(aln, coverage), as, "GRanges")
cvg <- plyr::ldply(cvg, as.data.frame)
colnames(cvg)[1] <- "iid"
cvg$panel <- cvg$iid

## discover splice junctions
sj <- lapply(aln, summarizeJunctions)

## ignore splices to out-of-range places
sj <- lapply(sj, subsetByOverlaps, zoom, type = "within")

## make swoop shapes representing splicing events
sj <- plyr::ldply(sj, as.data.frame)
swoops <- make.swoop(sj, scale = splice.scale)
swoops$panel <- swoops$.id

## add metadata, if any
if (!is.null(meta)){
  cvg <- merge(cvg, meta, all.x = TRUE)
}

## rescale by total coverage, if provided
if (!is.null(meta$reads)){
  cvg$depth <- with(cvg, score/(reads/1e6))
}else{
  cvg$depth <- cvg$score
}

## set colours for coverage panels, if not provided
if (is.null(colours)) {
  colours <- rep("darkgrey", length(aln))
  names(colours) <- names(aln)
}

## set x-axis limits to cover region of interest
xlims <- range(start(zoom), end(zoom))

## draw coverage plot with splice events
cvg <- subset(cvg, depth > min.coverage)
cvg$depth <- pmin(cvg$depth, max.coverage)
swoops <- subset(swoops, score > min.splices)

# whole gene
start_lim <- 96070333
end_lim <- 96161547

# # zoom
# start_lim <- 96069079
# end_lim <- 96110000

# plot
p0 <- 
  ggplot() +
  geom_rect(data = cvg, aes(xmin = start, xmax = end, ymin = 0, ymax = depth, fill = iid)) +
  geom_rect(data = del_gr, aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5), fill = "orange", alpha = 0.3) + 
  geom_line(data = swoops, aes(x = x, y = y, group = row, colour = .id)) +
  scale_x_continuous(limits = xlims) +
  scale_y_continuous(limits = c(min(swoops$y), 0.5)) +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours) +
  guides(fill = FALSE, colour = FALSE) +
  facet_grid(panel ~ ., scale = "free") +
  ylab("coverage (reads per million)\n") +
  coord_cartesian(xlim = c(start_lim, end_lim)) +
  theme_gbrowse() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())

# ggsave(filename = "test.pdf", plot = p0, width = 15, height = 10)

## draw transcript ideograms
p1 <-
  ggplot() +
  geom_tx(tx[tx$type == "exon"], at = 0, height = 5, fill = "grey50", colour.by = colour.by) +
  scale_x_continuous("\nposition (Mb)", limits = xlims, labels = function(x) x/1e6) +
  scale_fill_manual(values = c("grey50", "black", "blue")) +
  facet_grid(panel ~ .) +
  guides(fill = FALSE) +
  ylab("") +
  coord_cartesian(xlim = c(start_lim, end_lim)) +
  theme_gbrowse() +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

## combine plots, plot
pdf("sashimi_plot.CNOT6L.KO.all.pdf", width = 15, height = 10)
cowplot::plot_grid(p0, p1, nrow = 2, rel_heights = c(3/4, 1/4))
dev.off()


