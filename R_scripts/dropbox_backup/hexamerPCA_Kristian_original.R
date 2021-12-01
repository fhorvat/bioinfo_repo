library(Biostrings)
library(ggplot2)

setwd("C:/Users/kristian/Dropbox/Petr/MT")

ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
makeNewSample <- FALSE

if(! makeNewSample) {
  tsHist <- "20160801-223646"
}


comboClasses <- list(
  "MLT1" = c("MLT1A", "MLT1A0", "MLT1A1", "MLT1B", "MLT1C", "MLT1D", "MLT1E", "MLT1E1", "MLT1E1A",
             "MLT1E2", "MLT1E3", "MLT1F", "MLT1F1", "MLT1F2", "MLT1G", "MLT1G1", "MLT1G3", "MLT1H",
             "MLT1H1", "MLT1H2", "MLT1I", "MLT1J", "MLT1J1", "MLT1J2", "MLT1K", "MLT1L", "MLT1M",
             "MLT1N2", "MLT1O"),
  "MLT2" = c("MLT2B1", "MLT2B2", "MLT2B3", "MLT2B4", "MLT2B5", "MLT2C1", "MLT2C2", "MLT2D", "MLT2E",
             "MLT2F"),
#  "MLTR"  = c("MLTR11A", "MLTR11B", "MLTR12", "MLTR13", "MLTR14", "MLTR18_MM", "MLTR18A_MM",
#             "MLTR18B_MM", "MLTR18C_MM", "MLTR18D_MM", "MLTR25A", "MLTR25C", "MLTR31A_MM",
#             "MLTR31C_MM", "MLTR31D_MM", "MLTR31E_MM", "MLTR31F_MM", "MLTR31FA_MM", "MLTR32C_MM",
#             "MLTR73"),
  "MT2"   = c("MT2_Mm"),
  "MT2A"  = c("MT2A"),
  "MT2B"  = c("MT2B", "MT2B1", "MT2B2"),
  "MT2C"  = c("MT2C_Mm"),
  "MTA"   = c("MTA_Mm"),
  "MTB"   = c("MTB", "MTB_Mm"),
#  "MTB2"  = c("MTB_Mm"),
  "MTC"   = c("MTC"),
  "MTD"   = c("MTD"),
  "MTE"   = c("MTEa", "MTEb"),
  "MTE2"  = c("MTE2a", "MTE2b"),
  "ORR1A" = c("ORR1A0", "ORR1A1", "ORR1A2", "ORR1A3", "ORR1A4"),
  "ORR1B" = c("ORR1B1", "ORR1B2"),
  "ORR1C" = c("ORR1C1", "ORR1C2"),
  "ORR1D" = c("ORR1D1", "ORR1D2"),
  "ORR1E" = c("ORR1E"),
  "ORR1F" = c("ORR1F"),
  "ORR1G" = c("ORR1G")
)

transClasses <- sapply(names(comboClasses), function(x) {
  data.frame(from = comboClasses[[x]], to = rep.int(x, length(comboClasses[[x]])), stringsAsFactors = F)
  }, simplify = F)
transClasses <- do.call(rbind, transClasses)
rownames(transClasses) <- transClasses$from

mtFile <- "MT_LTR_filtered_by_length.fasta"
mtFile <- "MT_MT2_ORR_MLT_fullLengthLTR.fasta"
mtFile <- "MT_MT2_ORR_MLT_properFullLengthLTR.fasta"
mtFile <- "retrotransposon_consensusLength_5percent_LTR.fasta"
MTs <- readDNAStringSet(mtFile)
mtNames <- as.data.frame(do.call(rbind, strsplit(names(MTs), split = "\\|")))
names(mtNames) <- c("coord", "MT")
mtNames$finalClass <- transClasses[as.character(mtNames$MT), "to"]

classInstances <- 200

totalsInClass <- table(mtNames$finalClass)
overThreshold <- names(totalsInClass[totalsInClass > classInstances])
richMTs <- MTs[mtNames$finalClass %in% overThreshold]
richMTNames <- mtNames[mtNames$finalClass %in% overThreshold,]

nClasses <- length(overThreshold)

if(makeNewSample) {
  mtSubset <- unlist(tapply(1:length(richMTs), mtNames$finalClass[mtNames$finalClass %in% overThreshold], function(x) x[sample.int(length(x), classInstances, replace = F)]))

  sampledMTs <- richMTs[mtSubset,]
  writeXStringSet(sampledMTs,
                  file=paste("Sampled_MTs_all_classes_", ts, ".fasta", sep = ""),
                  compress=TRUE)
} else {
  sampledMTs <- readDNAStringSet(file=paste("Sampled_MTs_all_classes_", tsHist, ".fasta", sep = ""))
}

mtKmerFreqs      <- oligonucleotideFrequency(sampledMTs, 6, as.prob = TRUE, with.labels = TRUE)
mtKmerFreqsTotal <- oligonucleotideFrequency(sampledMTs, 6, as.prob = TRUE, simplify.as = "collapsed")

rownames(mtKmerFreqs) <- richMTNames$finalClass[mtSubset]

hexamers <- colnames(mtKmerFreqs)
hexamersDss <- DNAStringSet(hexamers)
rcHexamers <- as.character(reverseComplement(hexamersDss))
identicalKmers <- sapply(hexamers, function(x) which(x == rcHexamers))
uniqueWords <- names(identicalKmers[identicalKmers > 2048])

revComs <- data.frame(sense = uniqueWords, stringsAsFactors = F)
revComs$anti <- as.character(reverseComplement(DNAStringSet(revComs$sense)))


goPalindrome <- sapply(hexamersDss, findPalindromes, min.armlength = 3, min.looplength = 0, max.looplength = 0, max.mismatch = 0)
isPalindrome <- sapply(goPalindrome, function(x) length(start(x)) > 0)
names(isPalindrome) <- hexamers
palindromes <- names(isPalindrome[isPalindrome])


hmTmp <- hexamers
revComs <- data.frame()
while (length(hmTmp) > 0) {
  com <- hmTmp[1]
  rcom <- as.character(reverseComplement(DNAStringSet(com)))

  revComs <- rbind(revComs, data.frame(sense = com, anti = rcom))
  hmTmp <- hmTmp[! hmTmp %in% c(com, rcom)]

}

if(makeNewSample) {
  mtClean <- mtKmerFreqs[,revComs$sense] + mtKmerFreqs[,revComs$anti]
  #mtClean[, palindromes[1:32]] <- mtClean[, palindromes[1:32]]/2
  save(mtClean, file=paste("mtClean_", ts, ".Robj", sep = ""))
  
  mtCleanUnique <- unique(mtClean)
  save(mtCleanUnique, file=paste("mtCleanUnique_", ts, ".Robj", sep = ""))
  
} else {
  load(file=paste("mtClean_", tsHist, ".Robj", sep = ""))
  load(file=paste("mtCleanUnique_", tsHist, ".Robj", sep = ""))
}

#mtCleanUniqueClass <- as.data.frame(mtCleanUnique)
mtCleanUniqueClass <- as.data.frame(mtClean)
mtCleanUniqueClass$class <- as.factor(rownames(mtCleanUniqueClass))

# tweak to remove MTLR
# mtCleanUniqueClass <- mtCleanUniqueClass[mtCleanUniqueClass$class != "MTLR",]

# mtCleanSubset <- mtCleanUniqueClass[sample.int(nrow(mtCleanUniqueClass), 500),]
#
# mtCleanSubset <- by(mtCleanUniqueClass, mtCleanUniqueClass$class, function(x) x[sample.int(nrow(x), 100),])
# mtCleanSubset <- do.call(rbind, mtCleanSubset)

if(makeNewSample) {
  mtRF <- randomForest::randomForest(class ~ ., mtCleanUniqueClass, importance = T)
  save(mtRF, file=paste("mtRF_", ts, ".Robj", sep = ""))
} else {
  load(file=paste("mtRF_", tsHist, ".Robj", sep = ""))
}

library(pheatmap)
#pheatmap(mtRFunclass$proximity)
pheatmap(mtRF$confusion[1:nClasses,1:nClasses], cluster_rows = F, cluster_cols = F)
pheatmap(mtRF$confusion[1:nClasses,1:nClasses], cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white", "black"))(100))

heatmapOrder <- hclust(dist(mtRF$confusion[1:nClasses,1:nClasses]))$order

#heatmapOrderManual <- c("MT2", "MT2C", "MT2B", "MT2A", "MLTR", "MLT2", "MLT1", "ORR1G", "ORR1F",
#                        "ORR1E", "ORR1D", "MTE2", "MTE", "ORR1C", "ORR1B", "ORR1A", "MTD", "MTC",
#                        "MTB2", "MTB", "MTA")
#heatmapOrderManual <- 1:nClasses
#heatmapOrderManual <- c("MLT1", "MLT2", "MLTR", "MT2A", "MT2B", "MT2C", "MT2", "MTA", "MTB", "MTB2",
#                        "MTC", "MTD", "MTE", "MTE2", "ORR1A", "ORR1B", "ORR1C", "ORR1D", "ORR1E",
#                        "ORR1F", "ORR1G")
#heatmapOrderManual <- c("MLT1", "MLT2", "MLTR", "MT2A", "MT2B", "MT2C", "MT2", "MTA", "MTB",
#                        "MTC", "MTD", "MTE", "MTE2", "ORR1A", "ORR1B", "ORR1C", "ORR1D", "ORR1E",
#                        "ORR1F", "ORR1G")
heatmapOrderManual <- c("MLT1", "MLT2", "MT2A", "MT2B", "MT2C", "MT2", "MTA", "MTB",
                        "MTC", "MTD", "MTE", "MTE2", "ORR1A", "ORR1B", "ORR1C", "ORR1D", "ORR1E",
                        "ORR1F", "ORR1G")

wgb <- colorRampPalette(c("white", "gray10", "black"))
if(makeNewSample) pdf(file=paste("MT_confusionmatrix_", ts, ".pdf", sep = ""), width=9, height=8)
pheatmap(mtRF$confusion[heatmapOrderManual, heatmapOrderManual],
         treeheight_row = 10,
         treeheight_col = 10,
         cluster_rows = F, cluster_cols = F,
#         display_numbers = T,
         color = wgb(100))
if(makeNewSample) dev.off()

library(corrplot)
corrplot(mtRF$confusion[heatmapOrder, heatmapOrder],
         is.corr = F,
         method = "color",
         order = "AOE",
         tl.col = "black",
         addCoef.col = "red",
         cl.lim = c(0,200)
         )

library(DECIPHER)
goodGuys <- sampledMTs[which(mtRF$y == mtRF$predicted)]
goodGuys <- OrientNucleotides(goodGuys)
goodGuysClass <- mtRF$y[which(mtRF$y == mtRF$predicted)]

goodGuysSample <- tapply(goodGuys, goodGuysClass, function(x) x[sample.int(length(x), 20)])
goodGuysToAlign <- unlist(DNAStringSetList(unlist(goodGuysSample)))

if(makeNewSample) {
  alignment <- AlignSeqs(OrientNucleotides(goodGuysToAlign), iterations = 5, refinements = 5)
  writeXStringSet(alignment, file=paste("All_classes_alignment_", ts, ".fasta", sep = ""), compress=FALSE)
} else {
  alignment <- readDNAStringSet(file=paste("All_classes_alignment_", tsHist, ".fasta", sep = ""))
}


ts <- tsHist
makeNewSample <- TRUE
goodGuysSample <- split(goodGuys, goodGuysClass)

sapply(goodGuysSample, function(x) summary(width(x)))

perGroupAlignment <- sapply(names(goodGuysSample), function(x) {
  if(makeNewSample) {
    oriented <- OrientNucleotides(goodGuysSample[[x]])
    as <- AlignSeqs(oriented, iterations = 5, refinements = 5, gapOpening=c(-56, -52), gapExtension=c(-35, -30), gapPower = 2)
    names(as) <- paste(x, names(as), sep=".")
    writeXStringSet(as, file=paste(x, "_allgood_alignment_", ts, ".fasta", sep = ""), compress=FALSE)
  } else {
    as <- readDNAStringSet(file=paste(x, "_allgood_alignment_", tsHist, ".fasta", sep = ""))
  }
    as
  })

if(makeNewSample) {
  cumulativeAlignment <- AlignProfiles(perGroupAlignment[[1]], perGroupAlignment[[2]])
  for(i in 3:length(perGroupAlignment)) {
    cumulativeAlignment <- AlignProfiles(cumulativeAlignment, perGroupAlignment[[i]])
  }
  writeXStringSet(cumulativeAlignment, file=paste("All_classes_cumlative_stepwise_alignment_", ts, ".fasta", sep = ""), compress=FALSE)
} else {
  cumulativeAlignment <- readDNAStringSet(file=paste("All_classes_cumlative_stepwise_alignment_", tsHist, ".fasta", sep = ""))
}
  
JS <- function(x, y){
  # Function to compute Shannon-Jensen Divergence
  # x and y are the frequencies for the same p categories
  # Assumes relative abundance transformation already happened (for efficiency)

  # Define the mean point
  m <- (x+y)/2
  # Define each samples component
  P1 <- x*log(x/m)
  P2 <- y*log(y/m)
  # In the case of zeroes entries log is undefined, JSD is defined as zero
  P1[!is.finite(P1)] <- 0
  P2[!is.finite(P2)] <- 0
  d <- (P1+P2)/2
  return(sum(d, na.rm = TRUE))
}

conservationProfiles <- sapply(perGroupAlignment, function(x) {

  baseFreq <- colSums(letterFrequency(x, "ACGT", OR = 0))
  baseFreq <- baseFreq / sum(baseFreq)
  baseFreq <- c(.25, .25, .25, .25)

  JSmax <- JS(c(1, 0, 0, 0), c(.25, .25, .25, .25))

  mal <- DNAMultipleAlignment(x)
  cm <- consensusMatrix(mal)
  cMatrix <- cm[c("A", "C", "G", "T"),]
  gMatrix <- cm["-", ]

  totSums <- colSums(rbind(cMatrix, gMatrix))

  freqMatrix <- cMatrix/totSums
  freqGap    <- 1 - gMatrix/totSums

  JSD <- apply(freqMatrix, 2, function(x) JS(x, baseFreq))
  JSDcorr <- JSD * freqGap / JSmax

#  cbRes <- rbind(cMatrix, gMatrix, JSDcorr)
#  cbRes
  
#  list(pos = 1:length(JSDcorr), cor = JSDcorr)
  
  
  JSDcorr
  
})

library(reshape2)
moltenConservationProfiles <- melt(conservationProfiles, value.name = "rawcons", level = "TR")
moltenConservationProfiles$conservation <- unlist(
  tapply(moltenConservationProfiles$rawcons, moltenConservationProfiles$LTR, 
         function(x)  {
           n <- length(x) 
           runmed(x, k = 1 + 2 * min((n-1)%/% 2, ceiling(0.025*n)))
         }
  )
)
moltenConservationProfiles$seq <- unlist(tapply(moltenConservationProfiles$LTR, moltenConservationProfiles$LTR, function(x) 1:length(x)))

ggplot(moltenConservationProfiles, aes(x=seq, y=conservation, color=LTR, group=LTR)) +
  geom_line() +
  facet_wrap(~ LTR, ncol=2, scales = "free_x")
# +
#  geom_smooth(se = FALSE, span = 0.015)

#mtCleanUniqueNorm <- log(t(t(mtCleanUnique)/apply(mtCleanUnique,1, function(x)exp(mean(log(x+1))))))

hexVariances <- sapply(as.data.frame(mtCleanUnique), var)

mtTsne <- Rtsne(mtCleanUnique)

