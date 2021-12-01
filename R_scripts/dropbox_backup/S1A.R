rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd(".")

######################################################## LIBRARIES
library(Biostrings)
library(randomForest)
library(pheatmap)
library(ggplot2)

######################################################## FUNCTIONS

######################################################## READ DATA
comboClasses <- list(
  "MLT1" = c("MLT1A", "MLT1A0", "MLT1A1", "MLT1B", "MLT1C", "MLT1D", "MLT1E", "MLT1E1", "MLT1E1A",
             "MLT1E2", "MLT1E3", "MLT1F", "MLT1F1", "MLT1F2", "MLT1G", "MLT1G1", "MLT1G3", "MLT1H",
             "MLT1H1", "MLT1H2", "MLT1I", "MLT1J", "MLT1J1", "MLT1J2", "MLT1K", "MLT1L", "MLT1M",
             "MLT1N2", "MLT1O"),
  "MLT2" = c("MLT2B1", "MLT2B2", "MLT2B3", "MLT2B4", "MLT2B5", "MLT2C1", "MLT2C2", "MLT2D", "MLT2E",
             "MLT2F"),
  "MT2"   = c("MT2_Mm"),
  "MT2A"  = c("MT2A"),
  "MT2B"  = c("MT2B", "MT2B1", "MT2B2"),
  "MT2C"  = c("MT2C_Mm"),
  "MTA"   = c("MTA_Mm"),
  "MTB"   = c("MTB", "MTB_Mm"),
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

# sequences of LTRs with length +-5% of consensus length
mtFile <- "retrotransposon_consensusLength_5percent_LTR.fasta"
MTs <- readDNAStringSet(mtFile)
mtNames <- as.data.frame(do.call(rbind, strsplit(names(MTs), split = "\\|")))
names(mtNames) <- c("coord", "MT")
mtNames$finalClass <- transClasses[as.character(mtNames$MT), "to"]

######################################################## MAIN CODE
# get sequences of 200 random LTRs per class
classInstances <- 200
totalsInClass <- table(mtNames$finalClass)
overThreshold <- names(totalsInClass[totalsInClass > classInstances])
richMTs <- MTs[mtNames$finalClass %in% overThreshold]
richMTNames <- mtNames[mtNames$finalClass %in% overThreshold, ]
nClasses <- length(overThreshold)
mtSubset <- unlist(tapply(1:length(richMTs), mtNames$finalClass[mtNames$finalClass %in% overThreshold], function(x) x[sample.int(length(x), classInstances, replace = F)]))
sampledMTs <- richMTs[mtSubset, ]

# hexamers frequency in LTR sequences
mtKmerFreqs <- oligonucleotideFrequency(sampledMTs, 6, as.prob = TRUE, with.labels = TRUE)
mtKmerFreqsTotal <- oligonucleotideFrequency(sampledMTs, 6, as.prob = TRUE, simplify.as = "collapsed")
rownames(mtKmerFreqs) <- richMTNames$finalClass[mtSubset]

# unique hexamers
hexamers <- colnames(mtKmerFreqs)
hmTmp <- hexamers
revComs <- data.frame()
while(length(hmTmp) > 0) {
  com <- hmTmp[1]
  rcom <- as.character(reverseComplement(DNAStringSet(com)))
  revComs <- rbind(revComs, data.frame(sense = com, anti = rcom))
  hmTmp <- hmTmp[!hmTmp %in% c(com, rcom)]
}

# frequencies of unique hexameres in LTR sequences
mtClean <- mtKmerFreqs[, revComs$sense] + mtKmerFreqs[, revComs$anti]
mtCleanUniqueClass <- as.data.frame(mtClean)
mtCleanUniqueClass$class <- as.factor(rownames(mtCleanUniqueClass))

# random forest
mtRF <- randomForest::randomForest(class ~ ., mtCleanUniqueClass, importance = T)

# confusion matrix as heatmap
heatmapOrderManual <- c("MLT1", "MLT2", 
                        "ORR1G", "ORR1F", "ORR1E", "ORR1D", "ORR1C", "ORR1B", "ORR1A", 
                        "MTE", "MTE2", "MTD", "MTC", "MTB", "MTA", 
                        "MT2A", "MT2B", "MT2C", "MT2")
wgb <- colorRampPalette(c("white", "gray10", "black"))

pdf(file = "LTR_confusionMatrix.pdf", width = 9, height = 8, onefile = FALSE)
pheatmap(mtRF$confusion[heatmapOrderManual, heatmapOrderManual],
         treeheight_row = 10,
         treeheight_col = 10,
         cluster_rows = F, cluster_cols = F,
         color = wgb(100))
dev.off()