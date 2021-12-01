library(MUDAN)

## load built in 10X pbmcA dataset
data(pbmcA) 
pbmcA <- as.matrix(pbmcA) 
cd <- cleanCounts(pbmcA, min.reads = 10, min.detected = 10, verbose=FALSE)
mat <- normalizeCounts(cd, verbose=FALSE) 
matnorm.info <- normalizeVariance(mat, details=TRUE, verbose=FALSE) 
matnorm <- log10(matnorm.info$mat+1) 
pcs <- getPcs(matnorm[matnorm.info$ods,], nGenes=length(matnorm.info$ods), nPcs=30, verbose=FALSE) 
d <- dist(pcs, method='man')
emb <- Rtsne::Rtsne(d, is_distance=TRUE, perplexity=50, num_threads=parallel::detectCores(), verbose=FALSE)$Y 
rownames(emb) <- rownames(pcs)
com <- getComMembership(pcs, k=30, method=igraph::cluster_infomap, verbose=FALSE) 
dg <- getDifferentialGenes(cd, com)
plotEmbedding(emb, com, xlab=NA, ylab=NA, mark.clusters=TRUE, alpha=0.1, mark.cluster.cex=0.5, verbose=FALSE)


## Summarize gene expression within groups
## average expression
mat.summary <- do.call(cbind, lapply(levels(com), function(ct) {
  cells <- which(com==ct)
  if(length(cells) > 1) {
    Matrix::rowMeans(matnorm[, cells])
  } else {
    matnorm[,cells]
  }
}))
colnames(mat.summary) <- levels(com)
## fraction expressing
fe.summary <- do.call(cbind, lapply(levels(com), function(ct) {
  cells <- which(com==ct)
  if(length(cells) > 1) {
    Matrix::rowMeans(matnorm[, cells]>0)
  } else {
    matnorm[,cells]>0
  }
}))
colnames(fe.summary) <- levels(com)