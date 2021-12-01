library(GenomicAlignments)
flag0 <- scanBamFlag(isDuplicate = FALSE, isNotPassingQualityControls = FALSE)
param0 <- ScanBamParam(flag = flag0)

U3.GALP <- readGAlignmentPairs(samples_1cell_WT[1], use.names=TRUE, param = param0)

splice_overlaps <- findSpliceOverlaps(U3.GALP, exbytx)

head(U3.GALP)
head(first(U3.GALP))
head(last(U3.GALP))

U3.uqnames <- unique(names(U3.GALP))
U3.GALP_qnames <- factor(names(U3.GALP), levels=U3.uqnames)

U3.GALP_dup2unq <- match(U3.GALP_qnames, U3.GALP_qnames)

head(unique(cigar(first(U3.GALP))))

head(unique(cigar(last(U3.GALP))))

table(njunc(GenomicAlignments::first(U3.GALP)), njunc(GenomicAlignments::last(U3.GALP)))

colSums(cigarOpTable(cigar(GenomicAlignments::first(U3.GALP))))
colSums(cigarOpTable(cigar(GenomicAlignments::last(U3.GALP))))

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
TxDb.Mmusculus.UCSC.mm10.knownGene
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
exbytx <- exonsBy(txdb, by = "tx", use.names = TRUE)
length(exbytx) # nb of transcripts

exbytx_strand <- unlist(runValue(strand(exbytx)), use.names=FALSE)

tx <- transcripts(txdb, columns=c("tx_name", "gene_id"))
df <- mcols(tx)
exbytx2gene <- as.character(df$gene_id)
exbytx2gene <- factor(exbytx2gene, levels=unique(exbytx2gene))
names(exbytx2gene) <- df$tx_name
exbytx2gene <- exbytx2gene[names(exbytx)]
head(exbytx2gene)

nlevels(exbytx2gene)

U3.OV00 <- findOverlaps(U3.GALP, exbytx, ignore.strand=TRUE)
length(U3.OV00)

U3.GALP_ntx <- countQueryHits(U3.OV00)
mcols(U3.GALP)$ntx <- U3.GALP_ntx
head(U3.GALP)

table(U3.GALP_ntx)
mean(U3.GALP_ntx >= 1)

U3.GALP_ntx_again <- countOverlaps(U3.GALP, exbytx, ignore.strand=TRUE)
stopifnot(identical(unname(U3.GALP_ntx_again), U3.GALP_ntx))

U3.OV10 <- remapHits(U3.OV00, Lnodes.remapping=U3.GALP_qnames)
U3.uqnames_ntx <- countQueryHits(U3.OV10)
names(U3.uqnames_ntx) <- U3.uqnames
table(U3.uqnames_ntx)
mean(U3.uqnames_ntx >= 1)

U3.exbytx_nOV10 <- countSubjectHits(U3.OV10)
names(U3.exbytx_nOV10) <- names(exbytx)
mean(U3.exbytx_nOV10 >= 50)

head(sort(U3.exbytx_nOV10, decreasing=TRUE), n=10)

# encode
U3.grl <- grglist(U3.GALP)
U3.ovenc <- encodeOverlaps(U3.grl, exbytx, hits=U3.OV00, flip.query.if.wrong.strand=TRUE)
U3.ovenc

U3.unique_encodings <- levels(U3.ovenc)
length(U3.unique_encodings)
head(U3.unique_encodings)
U3.ovenc_table <- table(encoding(U3.ovenc))
tail(sort(U3.ovenc_table))

sort(U3.ovenc_table[isCompatibleWithSplicing(U3.unique_encodings)])
U3.OV00_is_comp <- isCompatibleWithSplicing(U3.ovenc)
table(U3.OV00_is_comp)

U3.compOV00 <- U3.OV00[U3.OV00_is_comp]
U3.GALP_ncomptx <- countQueryHits(U3.compOV00)
mcols(U3.GALP)$ncomptx <- U3.GALP_ncomptx

U3.noncompOV00 <- U3.OV00[!U3.OV00_is_comp]
U3.GALP_nnoncomptx <- countQueryHits(U3.noncompOV00)
mcols(U3.GALP)$nnoncomptx <- U3.GALP_nnoncomptx

# Number of compatible transcripts for each template:
U3.compOV10 <- remapHits(U3.compOV00, Lnodes.remapping=U3.GALP_qnames)
U3.uqnames_ncomptx <- countQueryHits(U3.compOV10)
names(U3.uqnames_ncomptx) <- U3.uqnames
table(U3.uqnames_ncomptx)
head(U3.GALP)

# Compute the reference query sequences
library(BSgenome.Mmusculus.UCSC.mm10)
U3.grl_first <- grglist(GenomicAlignments::first(U3.GALP, real.strand=TRUE), order.as.in.query=TRUE)
U3.grl_last <- grglist(GenomicAlignments::last(U3.GALP, real.strand=TRUE), order.as.in.query=TRUE)
U3.GALP_rqseq1 <- extractTranscriptSeqs(Mmusculus, U3.grl_first)
U3.GALP_rqseq2 <- extractTranscriptSeqs(Mmusculus, U3.grl_last)

U3.OV00_Lqstart <- extractQueryStartInTranscript(U3.grl, exbytx, hits=U3.OV00, ovenc=U3.ovenc)
head(subset(U3.OV00_Lqstart, U3.OV00_is_comp))

U3.OV00_Rqstart <- extractQueryStartInTranscript(U3.grl, exbytx, hits=U3.OV00, ovenc=U3.ovenc, for.query.right.end=TRUE)
head(subset(U3.OV00_Rqstart, U3.OV00_is_comp))

stepped_exons <- extractSteppedExonRanks(U3.ovenc) 
spanned_exons <- extractSpannedExonRanks(U3.ovenc)  
skipped_exons <- extractSkippedExonRanks(U3.ovenc) 

# splice overlaps 
splice_overlaps <- findSpliceOverlaps(U3.GALP, exbytx)
