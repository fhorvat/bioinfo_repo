### INFO: 
### DATE: 27. 11. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Deng_2014_Science_GSE45719/Analysis/velocyto/oocyte_no_2C")

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
# get color pallete
gg_color_hue <- function(n, times) {
  hues <- seq(15, 375, length = n + 1)
  hex_cl <- hcl(h = hues, l = 65, c = 100)[1:n]
  rep(hex_cl, times = times)
}

######################################################## READ DATA
# read sample table 
sample_df <- readr::read_delim(file = list.files(pattern = "*.sampleTable.tsv"), delim = "\t")

# named vector of colors
cell_colors <-
  sample_df$color %>%
  magrittr::set_names(., sample_df$sample_name)

# unique combination of colors and samples for legend
legend_df <- 
  sample_df %>% 
  dplyr::distinct(stage, color)

# read RDS
dat <- readRDS(file = list.files(pattern = "*.readsmartseq2dat.rds"))

######################################################## MAIN CODE
# ### Set up expression matrices, filtering genes to leave those that exceed some pre-defined g to the average expression magnitude
# # exonic read (spliced) expression matrix
# emat <- dat$emat
# 
# # intronic read (unspliced) expression matrix
# nmat <- dat$iomat
# 
# # spanning read (intron+exon) expression matrix
# smat <- dat$smat
# 
# # filter expression matrices based on some minimum max-cluster averages
# emat <- filter.genes.by.cluster.expression(emat, clusters = cell_colors, min.max.cluster.average = 5)
# nmat <- filter.genes.by.cluster.expression(nmat, clusters = cell_colors, min.max.cluster.average = 1)
# smat <- filter.genes.by.cluster.expression(smat, clusters = cell_colors, min.max.cluster.average = 0.5)
# 
# ### Several variants of velocity estimates using gene-relative model
# # kNN pooling with the gamma fit based on an extreme quantiles
# rvel.qf <- gene.relative.velocity.estimates(emat, nmat, kCells = 5, fit.quantile = 0.02, deltaT = 2)
# 
# # visualize
# pdf(file = "01.pca_velocity.gfit.kNN05.deltaT2.quantiles.pdf", width = 30, height = 10)
# pca.velocity.plot(rvel.qf, nPcs = 3, plot.cols = 3, cell.colors = ac(cell_colors, alpha = 1), cex = 3, pcount = 0.1, pc.multipliers = c(1, -1, 1))
# plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
# legend(x = "topleft", inset = 0, legend = legend_df$stage, col = ac(legend_df$color, alpha = 1), pch = 19, cex = 3, pt.cex = 5, horiz = F)
# dev.off()
# 
# # vizualize individual genes
# pdf(file = "test2.pdf", width = 40, height = 10)
# gene.relative.velocity.estimates(emat, nmat, deltaT = 1, kCells = 5, fit.quantile = 0.02, old.fit = rvel.qf, show.gene = 'Chga', cell.emb = emb, cell.colors = cell.colors)
# dev.off()
# 
# # Alternatively, we calculate gene-relative velocity, using k=5 cell kNN pooling, but now using entire range of expression to determine slope gamma,
# # and using spanning reads (smat) to fit the gene offsets.
# rvel <- gene.relative.velocity.estimates(emat, nmat, smat = smat, kCells = 5, deltaT = 2, min.nmat.emat.slope = 0.1, min.nmat.smat.correlation = 0.1)
# 
# # visualize
# pdf(file = "02.pca_velocity.gfit.kNN05.deltaT2.all.pdf", width = 30, height = 10)
# pca.velocity.plot(rvel, nPcs = 3, plot.cols = 3, cell.colors = ac(cell_colors, alpha = 1), cex = 3, pcount = 0.1, pc.multipliers = c(1, 1, 1))
# plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
# legend(x = "topleft", inset = 0, legend = legend_df$stage, col = ac(legend_df$color, alpha = 1), pch = 19, cex = 3, pt.cex = 5, horiz = F)
# dev.off()
# 
# # Here we calculate the most basic version of velocity estimates, using relative gamma fit, without cell kNN smoothing:
# rvel1 <- gene.relative.velocity.estimates(emat, nmat, deltaT = 1, deltaT2 = 1, kCells = 1)
# 
# # visualize
# pdf(file = "03.pca_velocity.knn0.relativeGfit.pdf", width = 30, height = 10)
# pca.velocity.plot(rvel1, nPcs = 3, plot.cols = 3, cell.colors = ac(cell_colors, alpha = 1), cex = 3, pcount = 0.1, pc.multipliers = c(1, -1, -1))
# plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
# legend(x = "topleft", inset = 0, legend = legend_df$stage, col = ac(legend_df$color, alpha = 1), pch = 19, cex = 3, pt.cex = 5, horiz = F)
# dev.off()


### Velocity estimate based on gene structure
## Genome-wide model fit:
# start with unfiltered matrices, as we can use more genes in these types of estimates
emat <- dat$emat
nmat <- dat$iomat
nmat["Zfp599", ] <- 0
smat <- dat$smat
emat <- filter.genes.by.cluster.expression(emat, cell_colors, min.max.cluster.average = 7)

## get global velocity estimates
gvel <- global.velcoity.estimates(emat, nmat, rvel, dat$base.df, smat = smat, 
                                  deltaT = 1, kCells = 5, kGenes = 15, kGenes.trim = 5, 
                                  min.gene.cells = 0, min.gene.conuts = 500)


## visualize in PCA space
pdf(file = "01.pca.velocity.pdf", width = 30, height = 10)
pca.velocity.plot(gvel, nPcs = 3, plot.cols = 3, cell.colors = ac(cell_colors, alpha = 1), cex = 3, pcount = 0.1, pc.multipliers = c(1, -1, -1))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(x = "topleft", inset = 0, legend = legend_df$stage, col = ac(legend_df$color, alpha = 1), pch = 19, cex = 3, pt.cex = 5, horiz = F)
dev.off()

## visualize in an tSNE space
# calculate tSNE clusters
perplex <- 30
pdf(NULL)
x <- tSNE.velocity.plot(rvel, nPcs = 15, cell.colors = ac(cell_colors, alpha = 1), cex = 3, 
                        perplexity = perplex, norm.nPcs = NA, pcount = 0.1, scale = 'log', do.par = F, 
                        arrow.scale = 3, arrow.lwd = 1)
dev.off()

# plot velocity in tSNE space
pdf(file = str_c("02.tsne.velocity.perplexity", perplex, ".pdf"), width = 20, height = 10)
par(mfrow = c(1, 2), mar = c(2.5, 2.5, 2.5, 1.5), mgp = c(2, 0.65, 0), cex = 0.85)
show.velocity.on.embedding.cor(x$current.emb, gvel, n = 100, scale = 'log', cell.colors = ac(cell_colors, alpha = 1), 
                               cex = 3, arrow.scale = 3, arrow.lwd = 1, do.par = F)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(x = "topleft", inset = 0, legend = legend_df$stage, col = ac(legend_df$color, alpha = 1), pch = 19, cex = 3, pt.cex = 5, horiz = F)
dev.off()

# plot velocity vector field:
pdf(file = str_c("03.tsne.velocity_field.perplexity", perplex, ".pdf"), width = 20, height = 10)
par(mfrow = c(1, 2), mar = c(2.5, 2.5, 2.5, 1.5), mgp = c(2, 0.65, 0), cex = 0.85)
show.velocity.on.embedding.cor(x$current.emb, gvel, n = 100, scale = 'log', cell.colors = ac(cell_colors, alpha = 1), 
                               cex = 3, arrow.scale = 3, arrow.lwd = 1, do.par = F, 
                               show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 20)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(x = "topleft", inset = 0, legend = legend_df$stage, col = ac(legend_df$color, alpha = 1), pch = 19, cex = 3, pt.cex = 5, horiz = F)
dev.off()

# plot cell trajectory modeling
pdf(file = str_c("04.tsne.cell_trajectory.perplexity", perplex, ".pdf"), width = 20, height = 10)
par(mfrow = c(1, 2), mar = c(2.5, 2.5, 2.5, 1.5), mgp = c(2, 0.65, 0), cex = 0.85)
show.velocity.on.embedding.eu(x$current.emb, gvel, n = 100, scale = 'log', cell.colors = ac(cell_colors, alpha = 1),
                              cex = 3, nPcs = 15, sigma = 2.5, show.trajectories = TRUE, diffusion.steps = 500,
                              n.trajectory.clusters = 15, ntop.trajectories = 3, embedding.knn = T, 
                              control.for.neighborhood.density = TRUE, n.cores = 20, do.par = F)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(x = "topleft", inset = 0, legend = legend_df$stage, col = ac(legend_df$color, alpha = 1), pch = 19, cex = 3, pt.cex = 5, horiz = F)
dev.off()

