library("DESeq2")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("geneplotter")
library("ggplot2")
library("TxDb.Mmusculus.UCSC.mm9.knownGene")

setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/")
filenames <- file.path(c("/common/WORK/vfranke/Projects/RNASeqEmbryo/Data/Mapped/Star_PairEnd/Star_PairEnd_EnsemblAnnot/s_1cell.PA/s_1cell.PA.bam",
                         "/common/WORK/vfranke/Projects/RNASeqEmbryo/Data/Mapped/Star_PairEnd/Star_PairEnd_EnsemblAnnot/s_MII.PA/s_MII.PA.bam", 
                         "/common/WORK/vfranke/Projects/RNASeqEmbryo/Data/Mapped/Star_PairEnd/Star_PairEnd_EnsemblAnnot/s_1cell.WE/s_1cell.WE.bam", 
                         "/common/WORK/vfranke/Projects/RNASeqEmbryo/Data/Mapped/Star_PairEnd/Star_PairEnd_EnsemblAnnot/s_GV.WE/s_GV.WE.bam", 
                         "/common/WORK/vfranke/Projects/RNASeqEmbryo/Data/Mapped/Star_PairEnd/Star_PairEnd_EnsemblAnnot/s_MII.WE/s_MII.WE.bam"))
filenames <- filenames[1:2]
bamfiles <- BamFileList(filenames, yieldSize = 2000000)
ebg <- exonsBy(TxDb.Mmusculus.UCSC.mm9.knownGene, by = "gene")
register(MulticoreParam())
se_Fugaku <- summarizeOverlaps(features = ebg, 
                               reads = bamfiles, 
                               mode = "Union", 
                               singleEnd = FALSE, 
                               ignore.strand = TRUE)

dds_Fugaku <- DESeqDataSet(se_Fugaku, design = ~1)
fpkm_Fugaku <- fpkm(dds_Fugaku, robust = F)
write.csv(fpkm_Fugaku, "fpkm_Fugaku.csv")

counts_Fugaku <- assay(se_Fugaku)
write.csv(counts_Fugaku, "counts_Fugaku.csv")

# writing DESeq2 differential expression data.frame-s
  # ratios:
  # 1C/MII PA
  # 1C/MII WE
  # MII/GV WE
  # samples: s_1cell.PA.bam s_MII.PA.bam s_1cell.WE.bam s_GV.WE.bam s_MII.WE.bam

dds_1C_MII_PA <- DESeqDataSet(se_Fugaku[, c("s_1cell.PA.bam", "s_MII.PA.bam")], design = ~1)
dds_1C_MII_PA <- DESeq(dds_1C_MII_PA)
write.csv(data.frame(results(dds_1C_MII_PA)), "DESeq2_1CvsMII_PA.csv")

dds_1C_MII_WE <- DESeqDataSet(se_Fugaku[, c("s_1cell.WE.bam", "s_MII.WE.bam")], design = ~1)
dds_1C_MII_WE <- DESeq(dds_1C_MII_WE)
write.csv(data.frame(results(dds_1C_MII_WE)), "DESeq2_1CvsMII_WE.csv")

dds_MII_GV_WE <- DESeqDataSet(se_Fugaku[, c("s_MII.WE.bam", "s_GV.WE.bam")], design = ~1)
dds_MII_GV_WE <- DESeq(dds_MII_GV_WE)
write.csv(data.frame(results(dds_MII_GV_WE)), "DESeq2_MIIvsGV_WE.csv")
