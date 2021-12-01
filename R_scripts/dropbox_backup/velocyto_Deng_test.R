### INFO:
### DATE: 27. 11. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Deng_2014_Science_GSE45719/Analysis/velocyto")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

library(velocyto.R)
library(XML)

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
gg_color_hue <- function(n, times) {
  hues <- seq(15, 375, length = n + 1)
  hex_cl <- hcl(h = hues, l = 65, c = 100)[1:n]
  rep(hex_cl, times = times)
}

######################################################## READ DATA
# list bam files
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Deng_2014_Science_GSE45719/Data/Mapped/STAR_mm10/4_bam_samples"
bam_files <-
  list.files(path = bam_path, pattern = "*.bam$", full.names = T, recursive = T) %>%
  .[!str_detect(., "liver|fibroblast")] %>% 
  magrittr::set_names(., str_replace(., ".*\\/(.*).SE.genome.Aligned.sortedByCoord.out.sample.bam", "\\1"))

### filter out smartSEQ2 samples
geo_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Deng_2014_Science_GSE45719/Data/Documentation/Deng_2014_Science_GSE45719.runInfo.xml"
geo_df <- 
  XML::xmlToDataFrame(geo_path, stringsAsFactors = F) %>% 
  as.tibble(.) %>% 
  dplyr::mutate(Ids = str_extract(Ids, "SRS\\d{6}")) %>% 
  dplyr::select(SRA_Sample_s = Ids, Description) %>% 
  dplyr::filter(str_detect(Description, "fibroblast|split|smartseq2")) %>% 
  dplyr::arrange(Description)

runinfo_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Deng_2014_Science_GSE45719/Data/Documentation/Deng_2014_Science_GSE45719.runInfo.txt"
runinfo_df <- 
  readr::read_delim(file = runinfo_path, delim = "\t") %>% 
  dplyr::select(Run_s, SRA_Sample_s) %>% 
  dplyr::filter(SRA_Sample_s %in% geo_df$SRA_Sample_s)

rename_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Deng_2014_Science_GSE45719/Data/Raw/Links/rename_links.sh"
rename_df <- 
  read_delim(file = rename_path, delim = " ", col_names = c("command", "run", "sample_name")) %>% 
  dplyr::mutate(Run_s = str_replace(run, ".fastq.gz", ""), 
                sample_name = str_replace(sample_name, ".SE.txt.gz", "")) %>% 
  dplyr::filter(Run_s %in% runinfo_df$Run_s)

# sample table with sample names and colors, order by stage
sample_df <- 
  tibble(sample_name = names(bam_files)) %>% 
  dplyr::filter(!(sample_name %in% rename_df$sample_name)) %>% 
  dplyr::mutate(stage = str_replace_all(sample_name, "_r.*|s_", ""),
                stage = str_replace(stage, "_", " "),
                color = gg_color_hue(n = length(unique(stage)), times = c(table(stage)))) %>% 
  dplyr::right_join(., 
                    tibble(stage = c("zygote", "early 2C_blastomere" , "mid 2C_blastomere", "late 2C_blastomere", 
                                     "2C blastomere", "4C blastomere", "8C blastomere", "16C blastomere",
                                     "early blastocyst", "mid blastocyst", "late blastocyst")), 
                    by = "stage")

# named vector of colors
cell_colors <- 
  sample_df$color %>% 
  magrittr::set_names(., sample_df$sample_name)

# unique combination of colors and samples for legend
legend_df <- 
  sample_df %>% 
  dplyr::distinct(stage, color)

# filter bam files
bam_files <- bam_files[names(bam_files) %in% sample_df$sample_name]

# refFlat annotation 
annot_path <- "/common/WORK/fhorvat/reference/mouse/mm10/UCSC/refFlat.mm10.20171127.txt.gz"

# parse gene annotation, annotate bam file reads
dat <- read.smartseq2.bams.modified(bam.files = bam_files, annotation.file = annot_path, n.cores = 10)

######################################################## MAIN CODE
### Set up expression matrices, filtering genes to leave those that exceed some pre-defined g to the average expression magnitude
# exonic read (spliced) expression matrix
emat <- dat$emat

# intronic read (unspliced) expression matrix
nmat <- dat$iomat

# spanning read (intron+exon) expression matrix
smat <- dat$smat

# filter expression matrices based on some minimum max-cluster averages
emat <- filter.genes.by.cluster.expression(emat, clusters = cell_colors, min.max.cluster.average = 5)
nmat <- filter.genes.by.cluster.expression(nmat, clusters = cell_colors, min.max.cluster.average = 1)
smat <- filter.genes.by.cluster.expression(smat, clusters = cell_colors, min.max.cluster.average = 0.5)

### Several variants of velocity estimates using gene-relative model
# kNN pooling with the gamma fit based on an extreme quantiles
rvel.qf <- gene.relative.velocity.estimates(emat, nmat, deltaT = 1, kCells = 5, fit.quantile = 0.02)

# visualize
pdf(file = "01_pca_velocity_kNN05_gfit_quantiles.pdf", width = 15, height = 10)
pca.velocity.plot(rvel.qf, nPcs = 6, plot.cols = 3, cell.colors = ac(cell_colors, alpha = 0.8), cex = 1.2, pcount = 0.1, pc.multipliers = c(1, -1, -1, -1, -1, 1))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(x = "topleft", inset = 0, legend = legend_df$stage, col = ac(legend_df$color, alpha = 0.8), lwd = 10, cex = 1, horiz = F)
dev.off()


# # vizualize individual genes
# pdf(file = "test2.pdf", width = 40, height = 10)
# gene.relative.velocity.estimates(emat, nmat, deltaT = 1, kCells = 5, fit.quantile = 0.02, old.fit = rvel.qf, show.gene = 'Chga', cell.emb = emb, cell.colors = cell.colors)
# dev.off()


# Alternatively, we calculate gene-relative velocity, using k=5 cell kNN pooling, but now using entire range of expression to determine slope gamma,
# and using spanning reads (smat) to fit the gene offsets.
rvel <- gene.relative.velocity.estimates(emat, nmat, smat = smat, deltaT = 1, kCells = 5, min.nmat.emat.slope = 0.1, min.nmat.smat.correlation = 0.1)

# visualize
pdf(file = "02_pca_velocity_kNN05_gfit_all.pdf", width = 15, height = 10)
pca.velocity.plot(rvel, nPcs = 6, plot.cols = 3, cell.colors = ac(cell_colors, alpha = 0.8), cex = 1.2, pcount = 0.1, pc.multipliers = c(1, -1, 1, 1, 1, 1))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(x = "topleft", inset = 0, legend = legend_df$stage, col = ac(legend_df$color, alpha = 0.8), lwd = 10, cex = 1, horiz = F)
dev.off()


# Here we calculate the most basic version of velocity estimates, using relative gamma fit, without cell kNN smoothing:
rvel1 <- gene.relative.velocity.estimates(emat, nmat, deltaT = 1, deltaT2 = 1, kCells = 1)

# visualize
pdf(file = "03_pca_velocity_gfit_relative.pdf", width = 15, height = 10)
pca.velocity.plot(rvel1, nPcs = 6, plot.cols = 3, cell.colors = ac(cell_colors, alpha = 0.8), cex = 1.2, pcount = 0.1, pc.multipliers = c(1, -1, 1, 1, 1, 1))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(x = "topleft", inset = 0, legend = legend_df$stage, col = ac(legend_df$color, alpha = 0.8), lwd = 10, cex = 1, horiz = F)
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
pdf(file = "04_pca_velocity_genome_wide.pdf", width = 15, height = 10)
pca.velocity.plot(gvel, nPcs = 6, plot.cols = 3, cell.colors = ac(cell_colors, alpha=0.8), cex = 1.2, pcount = 0.1, pc.multipliers = c(1, -1, -1, 1, 1, 1))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(x = "topleft", inset = 0, legend = legend_df$stage, col = ac(legend_df$color, alpha = 0.8), lwd = 10, cex = 1, horiz = F)
dev.off()

# visualize in tSNE space
perplex <- 20
pdf(file = str_c("05_tsne_velocity_genome_wide_perplex_", perplex, ".pdf"), width = 20, height = 10)
par(mfrow = c(1, 2), mar = c(2.5, 2.5, 2.5, 1.5), mgp = c(2, 0.65, 0), cex = 0.85)
x <- tSNE.velocity.plot(rvel, nPcs = 10, cell.colors = ac(cell_colors, alpha=0.8), cex = 1.5, perplexity = perplex, norm.nPcs = NA, pcount = 0.1, scale = 'log', do.par = F)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(x = "topleft", inset = 0, legend = legend_df$stage, col = ac(legend_df$color, alpha = 0.8), lwd = 10, cex = 1, horiz = F)
dev.off()


# Visualization on an existing embedding
# Here we use t-SNE embedding from the original publication (in emb variable).
vel <- rvel
emb <- x$projected.emb
pdf(file = "06_tsne_velocity_genome_wide_embeded.pdf", width = 20, height = 10)
par(mfrow = c(1, 2), mar = c(2.5, 2.5, 2.5, 1.5), mgp = c(2, 0.65, 0), cex = 0.85)
show.velocity.on.embedding.cor(emb, vel, n = 100, scale = 'sqrt', cell.colors = ac(cell_colors, alpha = 0.8), cex = 1.5, arrow.scale = 6, arrow.lwd = 1, do.par = F)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(x = "topleft", inset = 0, legend = legend_df$stage, col = ac(legend_df$color, alpha = 0.8), lwd = 10, cex = 1, horiz = F)
dev.off()

# Alternatively, the same function can be used to calculate a velocity vector field:
pdf(file = "07_tsne_velocity_vector_field.pdf", width = 20, height = 10)
par(mfrow = c(1, 2), mar = c(2.5, 2.5, 2.5, 1.5), mgp = c(2, 0.65, 0), cex = 0.85)
show.velocity.on.embedding.cor(emb, vel, n = 100, scale = 'sqrt', cell.colors = ac(cell_colors, alpha = 0.8), 
                               cex = 1, arrow.scale = 6, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 20, 
                               arrow.lwd = 1, do.par = F)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(x = "topleft", inset = 0, legend = legend_df$stage, col = ac(legend_df$color, alpha = 0.8), lwd = 10, cex = 1, horiz = F)
dev.off()

# Cell trajectory modeling
pdf(file = "08_cell_trajectory.pdf", width = 20, height = 10)
par(mfrow = c(1, 2), mar = c(2.5, 2.5, 2.5, 1.5), mgp = c(2, 0.65, 0), cex = 0.85)
x <- show.velocity.on.embedding.eu(emb, vel, n = 40, scale = 'sqrt', cell.colors = ac(cell_colors, alpha = 0.8),
                                   cex = 1, nPcs = 30, sigma = 2.5, show.trajectories = TRUE, diffusion.steps = 500,
                                   n.trajectory.clusters = 15, ntop.trajectories = 1, embedding.knn = T, 
                                   control.for.neighborhood.density = TRUE, n.cores = 20, do.par = F)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(x = "topleft", inset = 0, legend = legend_df$stage, col = ac(legend_df$color, alpha = 0.8), lwd = 10, cex = 1, horiz = F)
dev.off()
