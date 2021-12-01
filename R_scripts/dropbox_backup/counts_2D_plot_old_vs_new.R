library("DESeq2")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("geneplotter")
library("ggplot2")
library("TxDb.Mmusculus.UCSC.mm9.knownGene")

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

# filenames_old <- file.path(c("/common/WORK/vfranke/Projects/RNASeqEmbryo/Data/Mapped/Star_PairEnd/Star_PairEnd_EnsemblAnnot/s_1cell.PA/s_1cell.PA.bam",
#                              "/common/WORK/vfranke/Projects/RNASeqEmbryo/Data/Mapped/Star_PairEnd/Star_PairEnd_EnsemblAnnot/s_MII.PA/s_MII.PA.bam"))
# bamfiles <- BamFileList(filenames_old, yieldSize = 2000000)
# ebg <- exonsBy(TxDb.Mmusculus.UCSC.mm9.knownGene, by = "gene")
# register(MulticoreParam())
# se_PA <- summarizeOverlaps(features = ebg, 
#                            reads = bamfiles, 
#                            mode = "Union", 
#                            singleEnd = FALSE, 
#                            ignore.strand = TRUE)
# counts_old <- data.frame(assay(se_PA))

counts_old <- read.csv("counts_old.csv", row.names = 1)

counts_new <- byapply(counts_new, 3, rowMeans)

counts_all <- merge(counts_new, counts_old, by = 0, all = T)
rownames(counts_all) <- counts_all$Row.names
counts_all <- counts_all[, -1]
colnames(counts_all) <- c("GV_KO", "GV_WT", 
                          "MII_KO", "MII_WT", 
                          "oneC_KO", "oneC_WT", 
                          "oneC_PA", "MII_PA")

counts_all[counts_all == 0] <- NA
counts_all <- counts_all[rowSums(is.na(counts_all)) < 1, ]
counts_all_log <- log2(counts_all)

range_x <- range(counts_all_log$oneC_WT, counts_all_log$MII_WT)
range_y <- range(counts_all_log$oneC_PA, counts_all_log$MII_PA)

fit <- lm(oneC_WT ~ 0 + oneC_PA, data = counts_all_log)
lb <- paste("R^2 == ", format(summary(fit)$adj.r.squared, digits = 4))   
ggplot(counts_all_log, aes(x = oneC_WT, y = oneC_PA)) +
  geom_point(size = 1) +
  geom_abline(intercept = 0, slope = coef(fit)[1], color = "red") +
  annotate("text", x = range_x[2] - 1, y = range_y[1], label = lb, parse = TRUE) +
  scale_x_continuous("log2 counts 1C WT", limit = range_x) + 
  scale_y_continuous("log2 counts 1C PA", limit = range_y) +
  ggtitle("")
savepng("log2_counts_1C_WT_vs_1C_PA_2", width = 1000)

fit <- lm(MII_WT ~ 0 + MII_PA, data = counts_all_log)
lb <- paste("R^2 == ", format(summary(fit)$adj.r.squared, digits = 4))   
ggplot(counts_all_log, aes(x = MII_WT, y = MII_PA)) +
  geom_point(size = 1) +
  geom_abline(intercept = 0, slope = coef(fit)[1], color = "red") +
  annotate("text", x = range_x[2] - 1, y = range_y[1], label = lb, parse = TRUE) +
  scale_x_continuous("log2 counts MII WT", limit = range_x) + 
  scale_y_continuous("log2 counts MII PA", limit = range_y) +
  ggtitle("")
savepng("log2_counts_MII_WT_vs_MII_PA_2", width = 1000)


