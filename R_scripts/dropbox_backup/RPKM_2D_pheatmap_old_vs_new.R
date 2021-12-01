library("DESeq2")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("genefilter")
library("geneplotter")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("TxDb.Mmusculus.UCSC.mm9.knownGene")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

byapply <- function(x, by, fun, ...)
{
  # Create index list
  if (length(by) == 1)
  {
    nc <- ncol(x)
    split.index <- rep(1:ceiling(nc / by), each = by, length.out = nc)
  } else # 'by' is a vector of groups
  {
    nc <- length(by)
    split.index <- by
  }
  index.list <- split(seq(from = 1, to = nc), split.index)
  
  # Pass index list to fun using sapply() and return object
  sapply(index.list, function(i)
  {
    do.call(fun, list(x[, i], ...))
  })
}

counts_RPKM <- function(filenames, TxDB){
  bamfiles <- BamFileList(filenames, yieldSize = 2000000)
  ebg <- exonsBy(TxDB, by = "gene")
  register(MulticoreParam())
  se <- summarizeOverlaps(features = ebg, 
                          reads = bamfiles, 
                          mode = "Union", 
                          singleEnd = FALSE, 
                          ignore.strand = TRUE)
  
  # RPKM
  counts_se <- assay(se)
  lib_size <- colSums(counts_se) 
  ncounts <- t(t(counts_se) / lib_size)
  
  exon_lengths <- width(ebg) 
  exon_lengths_by_gene <- split(exon_lengths, names(exon_lengths))
  gene_lengths <- sapply(exon_lengths_by_gene, sum) 
  names(gene_lengths) <- sub("\\.\\d+", "", names(gene_lengths))
  
  common_genes <- intersect(row.names(ncounts), names(gene_lengths)) 
  subset_ncounts <- ncounts[row.names(ncounts) %in% common_genes, ] 
  gene_lengths <- gene_lengths[names(gene_lengths) %in% common_genes] 
  ncounts <- subset_ncounts/gene_lengths
  rpkm <- ncounts * 1e9 
  rpkm <- as.data.frame(rpkm)
  return(rpkm)  
}

# RPKM pheatmap function (no 0/NA/Inf, log2 values)
rpkm_pheatmap <- function(rpkm, main_title){
  rpkm_all_no_NA <- rpkm
  rpkm_all_no_NA[rpkm_all_no_NA == 0] <- NA
  rpkm_all_no_NA <- do.call(data.frame, lapply(rpkm_all_no_NA, function(x) replace(x, is.infinite(x), NA)))
  rownames(rpkm_all_no_NA) <- rownames(rpkm_all)
  rpkm_all_no_NA <- rpkm_all_no_NA[rowSums(is.na(rpkm_all_no_NA)) == 0, ]
  rpkm_all_no_NA <- log2(rpkm_all_no_NA)
  bk <- seq(range(rpkm_all_no_NA)[1], range(rpkm_all_no_NA)[2], length = 50)
  hmcols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(bk) - 1)
  pheatmap(rpkm_all_no_NA,
           clustering_method = "ward.D",
           main = main_title,
           col = hmcols, 
           breaks = bk,
           cluster_row = T, 
           cluster_cols = F,
           treeheight_row = 50, 
           show_rownames = F)
}

# new data RPKM - WT means all 3 timepoints
# a <- 1:18
# filenames_new <- file.path("/common/WORK/kristian/Projekti/Petr/Cnot6L/Mapping", paste0("11919X", a, "_sorted.bam"))
# rpkm_new <- counts_RPKM(filenames = filenames_new, TxDB = TxDb.Mmusculus.UCSC.mm10.knownGene) 
# colnames(rpkm_new) <- paste0(sampleTable$Time.Course[1:18], "_", sampleTable$Treatment.Control[1:18], rep(c(1, 2, 3), 6))
rpkm_new <- read.csv("rpkm_new.csv", row.names = 1)

# old data
# filenames_old <- file.path(c("/common/WORK/vfranke/Projects/RNASeqEmbryo/Data/Mapped/Star_PairEnd/Star_PairEnd_EnsemblAnnot/s_1cell.PA/s_1cell.PA.bam",
#                              "/common/WORK/vfranke/Projects/RNASeqEmbryo/Data/Mapped/Star_PairEnd/Star_PairEnd_EnsemblAnnot/s_MII.PA/s_MII.PA.bam"))
# rpkm_old <- counts_RPKM(filenames = filenames_old, TxDB = TxDb.Mmusculus.UCSC.mm9.knownGene) 
# colnames(rpkm_old) <- c("oneC_old", "MII_old")
rpkm_old <- read.csv("rpkm_old.csv", row.names = 1)

# new data WT only, timepoint means
rpkm_new_WT <- rpkm_new[, grep("WT", colnames(rpkm_new))]
rpkm_new_WT_means <- byapply(rpkm_WT, 3, rowMeans)
colnames(rpkm_new_WT_means) <- c("GV_triplicate_mean", "MII_triplicate_mean", "1C_triplicate_mean")

# pheatmap log RPKM WT data (new data: no timepoint means)
rpkm_all <- merge(rpkm_new_WT, rpkm_old, by = 0, all = T)
rownames(rpkm_all) <- rpkm_all$Row.names
rpkm_all <- rpkm_all[, c(-1:-4)]
colnames(rpkm_all) <- c("MII_WT1_new", "MII_WT2_new", "MII_WT3_new", "oneC_WT1_new", "oneC_WT2_new", "oneC_WT3_new", "oneC_old", "MII_old")
rpkm_all <- rpkm_all[, c(1, 2, 3, 8, 4, 5, 6, 7)]
rpkm_pheatmap(rpkm_all, main_title = "log2 RPKM WT old/new data")
savepng("pheatmap_log2_RPKM_old_new_WT_1", width = 1000)

# pheatmap log RPKM WT data (new data: timepoint means)
rpkm_all <- merge(rpkm_new_WT_means, rpkm_old, by = 0, all = T)
rownames(rpkm_all) <- rpkm_all$Row.names
rpkm_all <- rpkm_all[, c(-1, -2)]
colnames(rpkm_all) <- c("MII_new_triplicate_mean", "oneC_new_triplicate_mean", "oneC_old", "MII_old")
rpkm_all <- rpkm_all[, c(1, 4, 2, 3)]
rpkm_pheatmap(rpkm_all, main_title = "log2 RPKM WT old/new data")
savepng("pheatmap_log2_RPKM_old_new_WT_2", width = 1000)

# 2D plot log RPKM WT data (new data: timepoint means)
rpkm_all[rpkm_all == 0] <- NA
rpkm_all <- log2(rpkm_all)

ggplot(rpkm_all, aes(x = MII_new_triplicate_mean, y = MII_old)) +
  geom_point(size = 1) +
  scale_x_continuous("log2 RKPM MII WT new data (triplicate mean)") + 
  scale_y_continuous("log2 RPKM MII WT old data")
savepng("2D_log2_RPKM_old_new_MII_WT", width = 1000)

ggplot(rpkm_all, aes(x = oneC_new_triplicate_mean, y = oneC_old)) +
  geom_point(size = 1) +
  scale_x_continuous("log2 RKPM 1C WT new data (triplicate mean)") + 
  scale_y_continuous("log2 RPKM 1C WT old data")
savepng("2D_log2_RPKM_old_new_1C_WT", width = 1000)

