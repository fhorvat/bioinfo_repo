library("DESeq2")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")

rRNA_gtf_file <- "/common/WORK/fhorvat/reference/mm9/mm9_rRNA_20160412.gtf"
mm9_rRNA_TxDb <- makeTxDbFromGFF(rRNA_gtf_file, dataSource = "UCSC_track_browser", organism = "Mus musculus")
mm9_rRNA_ebg <- exonsBy(mm9_rRNA_TxDb, by = "tx")

filenames_Macfarlan <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/Other/MousePolyA/2cell/Macfarlan_2014/mapping/STAR", 
                                 paste0("SRR38562", 2:4, "_Aligned.sortedByCoord.out.bam"))
bamfiles <- BamFileList(filenames_Macfarlan, yieldSize = 2000000)
register(MulticoreParam())
se_Macfarlan <- summarizeOverlaps(features = mm9_rRNA_ebg, 
                                  reads = bamfiles, 
                                  mode = "Union", 
                                  singleEnd = TRUE, 
                                  ignore.strand = TRUE)
dds <- DESeqDataSet(se_Macfarlan, design = ~1)
fpkm_Macfarlan <- fpkm(dds)
colnames(fpkm_Macfarlan) <- paste0("SRR38562", 2:4)
fpkm_Macfarlan <- data.frame(fpkm_Macfarlan)
fpkm_Macfarlan$average <- rowMeans(fpkm_Macfarlan)

filenames_Hamazaki <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/Other/MousePolyA/2cell/Hamazaki_2015/mapping/STAR",
                                paste0("DRR02156", 3:6, "_Aligned.sortedByCoord.out.bam"))
bamfiles <- BamFileList(filenames_Hamazaki, yieldSize = 2000000)
register(MulticoreParam())
se_Hamazaki <- summarizeOverlaps(features = mm9_rRNA_ebg, 
                                 reads = bamfiles, 
                                 mode = "Union", 
                                 singleEnd = TRUE, 
                                 ignore.strand = TRUE)
dds <- DESeqDataSet(se_Hamazaki, design = ~1)
fpkm_Hamazaki <- fpkm(dds)
colnames(fpkm_Hamazaki) <- paste0("DRR02156", 3:6)
fpkm_Hamazaki <- data.frame(fpkm_Hamazaki)
fpkm_Hamazaki$average <- rowMeans(fpkm_Hamazaki)

