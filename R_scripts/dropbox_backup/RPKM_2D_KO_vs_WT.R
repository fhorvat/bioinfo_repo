library("ggplot2")
library("geneplotter")

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

rpkm_all <- read.csv("rpkm_new.csv", row.names = 1)

rpkm_GV <- rpkm_all[, grep("GV", colnames(rpkm_all))]
rpkm_GV <- byapply(rpkm_GV, 3, rowMeans)
colnames(rpkm_GV) <- c("GV_KO_mean", "GV_WT_mean")
rpkm_GV[rpkm_GV == 0] <- NA
rpkm_GV <- log2(rpkm_GV)
rpkm_GV <- as.data.frame(rpkm_GV)
ggplot(rpkm_GV, aes(x = GV_KO_mean, y = GV_WT_mean)) +
  geom_point(size = 1) +
  scale_x_continuous("log2 RKPM GV KO (triplicate mean)", limit = c(-15, 15)) + 
  scale_y_continuous("log2 RKPM GV WT (triplicate mean)", limit = c(-15, 15))
savepng("2D_log2_RPKM_GV_KO_vs_WT", width = 1000)

rpkm_MII <- rpkm_all[, grep("MII", colnames(rpkm_all))]
rpkm_MII <- byapply(rpkm_MII, 3, rowMeans)
colnames(rpkm_MII) <- c("MII_KO_mean", "MII_WT_mean")
rpkm_MII[rpkm_MII == 0] <- NA
rpkm_MII <- log2(rpkm_MII)
rpkm_MII <- as.data.frame(rpkm_MII)
ggplot(rpkm_MII, aes(x = MII_KO_mean, y = MII_WT_mean)) +
  geom_point(size = 1) +
  scale_x_continuous("log2 RKPM MII KO (triplicate mean)", limit = c(-15, 15)) + 
  scale_y_continuous("log2 RKPM MII WT (triplicate mean)", limit = c(-15, 15))
savepng("2D_log2_RPKM_MII_KO_vs_WT", width = 1000)

rpkm_oneC <- rpkm_all[, grep("1C", colnames(rpkm_all))]
rpkm_oneC <- byapply(rpkm_oneC, 3, rowMeans)
colnames(rpkm_oneC) <- c("oneC_KO_mean", "oneC_WT_mean")
rpkm_oneC[rpkm_oneC == 0] <- NA
rpkm_oneC <- log2(rpkm_oneC)
rpkm_oneC <- as.data.frame(rpkm_oneC)
ggplot(rpkm_oneC, aes(x = oneC_KO_mean, y = oneC_WT_mean)) +
  geom_point(size = 1) +
  scale_x_continuous("log2 RKPM 1C KO (triplicate mean)", limit = c(-15, 15)) + 
  scale_y_continuous("log2 RKPM 1C WT (triplicate mean)", limit = c(-15, 15))
savepng("2D_log2_RPKM_1C_KO_vs_WT", width = 1000)

