### INFO: 
### DATE: 27. 11. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/velocyto_files")

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

######################################################## READ DATA
# list bam files
files <- 
  list.files(path = outpath, pattern = "*unique.bam", full.names = T, recursive = T) %>% 
  magrittr::set_names(., str_replace(., ".*\\/(.*)_unique.bam", "\\1"))

# parse gene annotation, annotate bam file reads
dat <- read.smartseq2.bams(files, "genes.refFlat", n.cores = 10)

# read in cell cluster assignment and tSNE embedding used in the Furlan et al. (Science ’17).
cell.colors <- readRDS("cell.colors.rds")
emb <- readRDS("embedding.rds")

######################################################## MAIN CODE
### Gene filtering
# Spliced expression magnitude distribution across genes:
hist(log10(rowSums(dat$emat) + 1), col = 'wheat', xlab = 'log10[ number of reads + 1]', main = 'number of reads per gene')

### Set up expression matrices, filtering genes to leave those that exceed some pre-defined g to the average expression magnitude
# exonic read (spliced) expression matrix
emat <- dat$emat;

# intronic read (unspliced) expression matrix
nmat <- dat$iomat;

# spanning read (intron+exon) expression matrix
smat <- dat$smat;

# filter expression matrices based on some minimum max-cluster averages
emat <- filter.genes.by.cluster.expression(emat, cell.colors, min.max.cluster.average = 5)
nmat <- filter.genes.by.cluster.expression(nmat, cell.colors, min.max.cluster.average = 1)
smat <- filter.genes.by.cluster.expression(smat, cell.colors, min.max.cluster.average = 0.5)

# look at the resulting gene set
str(intersect(intersect(rownames(emat), rownames(nmat)), rownames(smat)))


### Several variants of velocity estimates using gene-relative model
# kNN pooling with the gamma fit based on an extreme quantiles
rvel.qf <- gene.relative.velocity.estimates(emat, nmat, deltaT = 1, kCells = 5, fit.quantile = 0.02)

# visualize
pdf(file = "test.pdf", width = 10, height = 10) 
pca.velocity.plot(rvel.qf, nPcs = 5, plot.cols = 2, cell.colors = ac(cell.colors, alpha = 0.7), cex = 1.2, pcount = 0.1, pc.multipliers = c(1, -1, -1, -1, -1))
dev.off()

# vizualize individual genes
pdf(file = "test2.pdf", width = 40, height = 10) 
gene.relative.velocity.estimates(emat, nmat, deltaT = 1, kCells = 5, fit.quantile = 0.02, old.fit = rvel.qf, show.gene = 'Chga', cell.emb = emb, cell.colors = cell.colors)
dev.off()


# Alternatively, we calculate gene-relative velocity, using k=5 cell kNN pooling, but now using entire range of expression to determine slope gamma, 
# and using spanning reads (smat) to fit the gene offsets.
rvel <- gene.relative.velocity.estimates(emat, nmat, smat = smat, deltaT = 1, kCells = 5, min.nmat.emat.slope = 0.1, min.nmat.smat.correlation = 0.1)

# visualize
pdf(file = "test3.pdf", width = 10, height = 10) 
pca.velocity.plot(rvel, nPcs = 5, plot.cols = 2, cell.colors = ac(cell.colors, alpha = 0.7), cex = 1.2, pcount = 0.1, pc.multipliers = c(1, -1, 1, 1, 1))
dev.off()


# Here we calculate the most basic version of velocity estimates, using relative gamma fit, without cell kNN smoothing:
rvel1 <- gene.relative.velocity.estimates(emat, nmat, deltaT = 1, deltaT2 = 1, kCells = 1)

# visualize
pdf(file = "test4.pdf", width = 10, height = 10) 
pca.velocity.plot(rvel1, nPcs = 5, plot.cols = 2, cell.colors = ac(cell.colors, alpha = 0.7), cex = 1.2, pcount = 0.1, pc.multipliers = c(1, -1, 1, 1, 1))
dev.off()


### Velocity estimate based on gene structure
## Genome-wide model fit:
# start with unfiltered matrices, as we can use more genes in these types of estimates
emat <- dat$emat
nmat <- dat$iomat
smat <- dat$smat
emat <- filter.genes.by.cluster.expression(emat, cell.colors, min.max.cluster.average = 7)
gvel <- global.velcoity.estimates(emat, nmat, rvel, dat$base.df, smat = smat, deltaT = 1, kCells = 5, kGenes = 15, kGenes.trim = 5, min.gene.cells = 0, min.gene.conuts = 500)

# visualize in PCA space
pdf(file = "test5.pdf", width = 10, height = 10) 
pca.velocity.plot(gvel, nPcs = 5, plot.cols = 2, cell.colors = ac(cell.colors, alpha=0.7), cex = 1.2, pcount = 0.1, pc.multipliers = c(1, -1, -1, 1, 1))
dev.off()

# visualize in tSNE space
pdf(file = "test6.pdf", width = 12, height = 6) 
par(mfrow = c(1, 2), mar = c(2.5, 2.5, 2.5, 1.5), mgp = c(2, 0.65, 0), cex = 0.85)
x <- tSNE.velocity.plot(rvel, nPcs = 15, cell.colors = cell.colors, cex = 0.9, perplexity = 200, norm.nPcs = NA, pcount = 0.1, scale = 'log', do.par = F)
dev.off()


# Visualization on an existing embedding
# Here we use t-SNE embedding from the original publication (in emb variable).
vel <- rvel
pdf(file = "test7.pdf", width = 9, height = 8) 
show.velocity.on.embedding.cor(emb, vel, n = 100, scale = 'sqrt', cell.colors = ac(cell.colors, alpha = 0.4), cex = 1, arrow.scale = 6, arrow.lwd = 1)
dev.off()

# Alternatively, the same function can be used to calculate a velocity vector field:
pdf(file = "test8.pdf", width = 9, height = 8) 
show.velocity.on.embedding.cor(emb, vel, n = 100, scale = 'sqrt', cell.colors = ac(cell.colors, alpha = 0.4), cex = 1, arrow.scale = 6, 
                               show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 20, arrow.lwd = 2)
dev.off()


# Cell trajectory modeling
pdf(file = "test9.pdf", width = 9, height = 8) 
x <- show.velocity.on.embedding.eu(emb, vel, n = 40, scale = 'sqrt', cell.colors = ac(cell.colors, alpha = 0.4),
                                   cex = 1, nPcs = 30, sigma = 2.5, show.trajectories = TRUE, diffusion.steps = 500,
                                   n.trajectory.clusters = 15, ntop.trajectories = 1, embedding.knn = T, control.for.neighborhood.density = TRUE, n.cores = 40)
dev.off()
