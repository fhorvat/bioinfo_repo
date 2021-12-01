### INFO: 
### DATE: 27. 11. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/velocyto/")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

library(velocyto.R)

######################################################## PATH VARIABLES
outpath <- getwd()

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))

######################################################## FUNCTIONS
# modified read.smartseq2.bams (in parallel::mclapply mc.preschedule = F)
read.smartseq2.bams.modified <- function (bam.files, annotation.file, min.exon.count = 100, n.cores = defaultNCores()){
  
  cat("reading gene annotation ... ")
  x <- read.delim(annotation.file, header = F, sep = "\t", stringsAsFactors = F)
  genes <- data.frame(name = x[, 1], chr = x[, 3], strand = x[, 4], start = x[, 5], end = x[, 6], stringsAsFactors = F)
  genes$p5 <- genes$start
  genes$p5[genes$strand == "-"] <- genes$end[genes$strand == "-"]
  genes$p3 <- genes$end
  genes$p3[genes$strand == "-"] <- genes$start[genes$strand == "-"]
  genes$size <- genes$end - genes$start
  cat("done (", nrow(genes), "genes)\n")
  
  cat("parsing exon information ... ")
  exons <- do.call(rbind, lapply(1:nrow(x), function(i) {
    df <- do.call(cbind, strsplit(as.character(x[i, c(10, 11)]), ","))
    cbind(df, rep(as.character(x[i, 1]), nrow(df)))
  }))
  exons <- data.frame(gene = as.character(exons[, 3]), start = as.numeric(exons[, 1]), end = as.numeric(exons[, 2]), stringsAsFactors = F)
  exons$chr <- genes$chr[match(exons$gene, genes$name)]
  exons <- exons[!duplicated(paste(exons[, 1], exons[, 2], exons[, 3])), ]
  genes <- genes[order(genes$size, decreasing = T), ]
  genes <- genes[!duplicated(genes$name), ]
  cat("done\n")
  
  cat("reading in", length(bam.files), "bam files ... ")
  cdl <- parallel::mclapply(bam.files, velocyto.R:::t.annotate.bam.reads,
                            genes = genes, exons = exons, margin = 1, exon.margin = 1,
                            mc.cores = n.cores, mc.preschedule = F)
  cat("done\n")
  
  cat("estimating gene counts ... ")
  edl <- parallel::mclapply(cdl, velocyto.R:::t.get.estimates2, genes = genes, mc.cores = n.cores)
  cat("done\n")
  
  cat("adjusting gene annotation based on expressed regions ... ")
  tl <- Matrix::colSums(do.call(rbind, parallel::mclapply(cdl, function(x) {
                                                            ect <- table(c(x$exonstart, x$exonend))
                                                            fect <- rep(0, nrow(exons))
                                                            fect[as.integer(names(ect))] <- ect
                                                            fect
                                                            }, mc.cores = n.cores, mc.preschedule = T)))
  expr.exons <- exons[tl > min.exon.count, ]
  expr.lstat <- velocyto.R:::lengthstats2(1e+08, genes = genes, exons = expr.exons)
  df <- data.frame(il = log10(expr.lstat[, 2] + 1), el = log10(expr.lstat[, 3] + 1))
  rownames(df) <- rownames(expr.lstat)
  df$nex <- as.integer(table(expr.exons$gene)[rownames(df)])
  df$nex[is.na(df$nex)] <- 0
  cat("done\n")
  
  emat <- do.call(cbind, lapply(edl, function(d) {
    (d[, "exon"])
  }))
  emat[!is.finite(emat)] <- 0
  emat <- as(emat, "dgCMatrix")
  smat <- do.call(cbind, lapply(edl, function(d) {
    (d[, "span"])
  }))
  smat[!is.finite(smat)] <- 0
  smat <- as(smat, "dgCMatrix")
  iomat <- do.call(cbind, lapply(edl, function(d) {
    (d[, "introno"])
  }))
  iomat[!is.finite(iomat)] <- 0
  iomat <- as(iomat, "dgCMatrix")
  return(list(emat = emat, iomat = iomat, smat = smat, base.df = df,
              exons = exons, genes = genes, expr.lstat = expr.lstat))
  
}

# get color pallete
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

######################################################## READ DATA
# list bam files
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10/velocyto_test"
bam_files <- 
  list.files(path = bam_path, pattern = "*.bam$", full.names = T, recursive = T) %>% 
  magrittr::set_names(., str_replace(., ".*\\/(.*).sample.bam", "\\1")) %>% 
  .[!str_detect(., "KO|Hamster")] 

# refFlat annotation 
annot_path <- "/common/WORK/fhorvat/reference/mouse/mm10/UCSC/refFlat.mm10.20171127.txt.gz"

# parse gene annotation, annotate bam file reads
dat <- read.smartseq2.bams.modified(bam.files = bam_files, annotation.file = annot_path, n.cores = 10)

# read in cell cluster assignment and tSNE embedding used in the Furlan et al. (Science ’17).
cell_colors <- readRDS("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/velocyto_files/cell.colors.rds")
emb <- readRDS("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/velocyto_files/embedding.rds")

######################################################## MAIN CODE
### set colors
# get unique color for each developmental stage
unique_cell_colors <- 
  tibble(sample = names(bam_files)) %>% 
  dplyr::mutate(sample = str_replace(sample, "_r.*", "")) %>% 
  dplyr::distinct(.) %>% 
  dplyr::mutate(color = gg_color_hue(nrow(.)))

# assign dev. stage color to each sample 
cell_colors_df <- 
  tibble(file_names = names(bam_files)) %>% 
  dplyr::mutate(sample = str_replace(file_names, "_r.*", "")) %>% 
  dplyr::left_join(., unique_cell_colors, by = "sample")

# create named vector of colors
cell_colors <- 
  cell_colors_df$color %>% 
  magrittr::set_names(., cell_colors_df$file_names)

  
### Gene filtering
# Spliced expression magnitude distribution across genes:
pdf(file = "01.reads_number.hist.pdf", width = 10, height = 10)
hist(log10(rowSums(dat$emat) + 1), col = 'wheat', xlab = 'log10[ number of reads + 1]', main = 'number of reads per gene')
dev.off()

### Set up expression matrices, filtering genes to leave those that exceed some pre-defined g to the average expression magnitude
# exonic read (spliced) expression matrix
emat <- dat$emat

# intronic read (unspliced) expression matrix
nmat <- dat$iomat

# spanning read (intron+exon) expression matrix
smat <- dat$smat

# filter expression matrices based on some minimum max-cluster averages
emat <- filter.genes.by.cluster.expression(emat, cell_colors, min.max.cluster.average = 5)
nmat <- filter.genes.by.cluster.expression(nmat, cell_colors, min.max.cluster.average = 1)
smat <- filter.genes.by.cluster.expression(smat, cell_colors, min.max.cluster.average = 0.5)

# look at the resulting gene set
str(intersect(intersect(rownames(emat), rownames(nmat)), rownames(smat)))


### Several variants of velocity estimates using gene-relative model
# kNN pooling with the gamma fit based on an extreme quantiles
rvel.qf <- gene.relative.velocity.estimates(emat, nmat, deltaT = 1, kCells = 5, fit.quantile = 0.02)

# visualize
pdf(file = "test.pdf", width = 10, height = 10) 
pca.velocity.plot(rvel.qf, nPcs = 5, plot.cols = 2, cell.colors = ac(cell_colors, alpha = 0.7), cex = 1.2, pcount = 0.1, pc.multipliers = c(1, -1, -1, -1, -1))
dev.off()

# vizualize individual genes
pdf(file = "test2.pdf", width = 10, height = 10) 
gene.relative.velocity.estimates(emat, nmat, deltaT = 1, kCells = 5, fit.quantile = 0.02, old.fit = rvel.qf, show.gene = "Ctnna3", cell.emb = NULL, cell.colors = cell_colors)
dev.off()


# Alternatively, we calculate gene-relative velocity, using k=5 cell kNN pooling, but now using entire range of expression to determine slope gamma, 
# and using spanning reads (smat) to fit the gene offsets.
rvel <- gene.relative.velocity.estimates(emat, nmat, smat = smat, deltaT = 1, kCells = 5, min.nmat.emat.slope = 0.1, min.nmat.smat.correlation = 0.1)

# visualize
pdf(file = "test3.pdf", width = 10, height = 10) 
pca.velocity.plot(rvel, nPcs = 5, plot.cols = 2, cell.colors = ac(cell_colors, alpha = 0.7), cex = 1.2, pcount = 0.1, pc.multipliers = c(1, -1, 1, 1, 1))
dev.off()


# Here we calculate the most basic version of velocity estimates, using relative gamma fit, without cell kNN smoothing:
rvel1 <- gene.relative.velocity.estimates(emat, nmat, deltaT = 1, deltaT2 = 1, kCells = 1)

# visualize
pdf(file = "test4.pdf", width = 10, height = 10) 
pca.velocity.plot(rvel1, nPcs = 5, plot.cols = 2, cell.colors = ac(cell_colors, alpha = 0.7), cex = 1.2, pcount = 0.1, pc.multipliers = c(1, -1, 1, 1, 1))
dev.off()


### Velocity estimate based on gene structure
## Genome-wide model fit:
# start with unfiltered matrices, as we can use more genes in these types of estimates
emat <- dat$emat
nmat <- dat$iomat
smat <- dat$smat
emat <- filter.genes.by.cluster.expression(emat, cell_colors, min.max.cluster.average = 7)
gvel <- global.velcoity.estimates(emat, nmat, rvel, dat$base.df, smat = smat, deltaT = 1, kCells = 5, kGenes = 15, kGenes.trim = 5, min.gene.cells = 0, min.gene.conuts = 500)

# visualize in PCA space
pdf(file = "test5.pdf", width = 10, height = 10) 
pca.velocity.plot(gvel, nPcs = 5, plot.cols = 2, cell.colors = ac(cell_colors, alpha=0.7), cex = 1.2, pcount = 0.1, pc.multipliers = c(1, -1, -1, 1, 1))
dev.off()

# visualize in tSNE space
pdf(file = "test6.pdf", width = 10, height = 10) 
par(mfrow = c(1, 1), mar = c(2.5, 2.5, 2.5, 1.5), mgp = c(2, 0.65, 0), cex = 0.85)
x <- tSNE.velocity.plot(rvel, 
                        nPcs = 3, 
                        cell.colors = cell_colors, 
                        cex = 0.9, 
                        perplexity = 2, 
                        norm.nPcs = NA, 
                        pcount = 0.1, 
                        scale = 'log', 
                        do.par = F)
dev.off()


# Visualization on an existing embedding
# Here we use t-SNE embedding from the original publication (in emb variable).
vel <- rvel
emb <- x$projected.emb
pdf(file = "test7.pdf", width = 10, height = 10)
show.velocity.on.embedding.cor(emb, vel, n = 100, scale = 'sqrt', cell.colors = ac(cell_colors, alpha = 0.4), cex = 1, arrow.scale = 6, arrow.lwd = 1)
dev.off()

# Alternatively, the same function can be used to calculate a velocity vector field:
pdf(file = "test8.pdf", width = 10, height = 10) 
show.velocity.on.embedding.cor(emb, vel, n = 5, scale = 'sqrt', cell.colors = ac(cell_colors, alpha = 0.4), cex = 1, arrow.scale = 6, 
                               show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 20, arrow.lwd = 2)
dev.off()


# Cell trajectory modeling
pdf(file = "test9.pdf", width = 10, height = 10) 
x <- show.velocity.on.embedding.eu(emb, vel, n = 5, scale = 'sqrt', cell.colors = ac(cell_colors, alpha = 0.4),
                                   cex = 1, nPcs = 3, sigma = 10, show.trajectories = TRUE, diffusion.steps = 500,
                                   n.trajectory.clusters = 3, ntop.trajectories = 3, embedding.knn = T, control.for.neighborhood.density = TRUE, n.cores = 40, 
                                   show.all.trajectories = T)
dev.off()

velocyto.R::pca.velocity.plot()
