library("DESeq2")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("TxDb.Mmusculus.UCSC.mm9.knownGene")

sampleTable <- read.csv("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/CNOT6L_sample_list_11919R_2015-10-29.csv", header = T)

TxDB <- TxDb.Mmusculus.UCSC.mm9.knownGene
a <- 1:18
filenames <- file.path("/common/WORK/kristian/Projekti/Petr/Cnot6L/Mapping/bbmap/mm9", paste0("11919X", a, "_sorted.bam"))
bamfiles <- BamFileList(filenames, yieldSize = 2000000)
ebg <- exonsBy(TxDB, by = "gene")

register(MulticoreParam())
se <- summarizeOverlaps(features = ebg, 
                        reads = bamfiles, 
                        mode = "Union", 
                        singleEnd = FALSE, 
                        ignore.strand = TRUE)

colData(se) <- DataFrame(sampleTable[1:18, ])
dds <- DESeqDataSet(se, design = ~1)
fpkm_CNOT6L <- fpkm(dds)
colnames(fpkm_CNOT6L) <- sampleTable$ID[1:18]
write.csv(fpkm_CNOT6L, "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/CNOT6L_mm9_fpkm.csv")
