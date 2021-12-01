library(ggbio)
library(GenomicRanges)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(biovizBase)
library(BSgenome.Hsapiens.UCSC.hg19)
library(VariantAnnotation)

p.ideo <- Ideogram(genome = "mm10")

## special highlights instead of zoomin!
png(filename = file.path(outpath, "ggbio.test2.png"), width = 1500, height = 400)
p.ideo + xlim(GRanges("chr15", IRanges(73131042, 73137993)))
dev.off()

# annotations
class(Homo.sapiens)
data(genesymbol, package = "biovizBase")
wh <- genesymbol[c("BRCA1", "NBR1")]
wh <- range(wh, ignore.strand = TRUE)

png(filename = file.path(outpath, "ggbio.test3.png"), width = 1500, height = 400)
autoplot(Homo.sapiens, which = wh)
dev.off()

png(filename = file.path(outpath, "ggbio.test4.png"), width = 1500, height = 400)
autoplot(Homo.sapiens, which = wh, label.color = "black", color = "brown", fill = "brown")
dev.off()

png(filename = file.path(outpath, "ggbio.test5.png"), width = 1500, height = 400)
autoplot(Homo.sapiens, which = wh, gap.geom = "chevron")
dev.off()

png(filename = file.path(outpath, "ggbio.test6.png"), width = 1500, height = 400)
autoplot(Homo.sapiens, which = wh, stat = "reduce")
dev.off()

png(filename = file.path(outpath, "ggbio.test7.png"), width = 1500, height = 400)
autoplot(Homo.sapiens, which = wh, columns = c("TXNAME", "GO"), names.expr = "TXNAME::GO")
dev.off()

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
png(filename = file.path(outpath, "ggbio.test8.png"), width = 1500, height = 400)
autoplot(txdb, which = wh)
dev.off()

gr.txdb <- crunch(txdb, which = wh)
colnames(values(gr.txdb))[4] <- "model"
grl <- split(gr.txdb, gr.txdb$tx_id)
names(grl) <- sample(LETTERS, size = length(grl), replace = TRUE)
png(filename = file.path(outpath, "ggbio.test9.png"), width = 1500, height = 400)
autoplot(grl, aes(type = model))
ggplot() + geom_alignment(grl, type = "model")
dev.off()

# sequence
bg <- BSgenome.Hsapiens.UCSC.hg19
p.bg <- autoplot(bg, which = wh)
png(filename = file.path(outpath, "ggbio.test10.png"), width = 1500, height = 400)
p.bg
dev.off()

png(filename = file.path(outpath, "ggbio.test11.png"), width = 1500, height = 200)
p.bg + zoom(1/100)
dev.off()

png(filename = file.path(outpath, "ggbio.test12.png"), width = 1500, height = 200)
p.bg + zoom(1/1000)
dev.off()

png(filename = file.path(outpath, "ggbio.test13.png"), width = 1500, height = 200)
p.bg + zoom(1/2500)
dev.off()

png(filename = file.path(outpath, "ggbio.test14.png"), width = 1500, height = 200)
autoplot(bg, which = resize(wh, width = width(wh)/2000), geom = "segment")
dev.off()


# bam files
fl.bam <- system.file("extdata", "wg-brca1.sorted.bam", package = "biovizBase")
wh <- keepSeqlevels(wh, "chr17")
png(filename = file.path(outpath, "ggbio.test15.png"), width = 1500, height = 200)
autoplot(fl.bam, which = wh)
dev.off()

png(filename = file.path(outpath, "ggbio.test16.png"), width = 1500, height = 200)
autoplot(fl.bam, which = resize(wh, width = width(wh)/10), geom = "gapped.pair")
dev.off()

bg <- BSgenome.Hsapiens.UCSC.hg19
p.mis <- autoplot(fl.bam, bsgenome = bg, which = wh, stat = "mismatch")
png(filename = file.path(outpath, "ggbio.test17.png"), width = 1500, height = 200)
p.mis
dev.off()

png(filename = file.path(outpath, "ggbio.test18.png"), width = 1500, height = 200)
autoplot(fl.bam, method = "estimate")
dev.off()

png(filename = file.path(outpath, "ggbio.test19.png"), width = 1500, height = 200)
autoplot(fl.bam, method = "estimate", which = paste0("chr", 17:18), aes(y = log(..coverage..)))
dev.off()

### variant annotation
fl.vcf <- system.file("extdata", "17-1409-CEU-brca1.vcf.bgz", package="biovizBase")
vcf <- readVcf(fl.vcf, "hg19")
vr <- as(vcf[, 1:3], "VRanges")
vr <- renameSeqlevels(vr, value = c("17" = "chr17"))
gr17 <- GRanges("chr17", IRanges(41234400, 41234530))
p.vr <- autoplot(vr, which = wh)

png(filename = file.path(outpath, "ggbio.test20.png"), width = 1500, height = 200)
p.vr
dev.off()

png(filename = file.path(outpath, "ggbio.test21.png"), width = 1500, height = 200)
p.vr + xlim(gr17)
dev.off()

png(filename = file.path(outpath, "ggbio.test22.png"), width = 1500, height = 200)
p.vr + xlim(gr17) + zoom()
dev.off()

### Building your tracks
## tks <- tracks(p.ideo, mismatch = p.mis, dbSNP = p.vr, ref = p.bs, gene = p.txdb)
## tks <- tracks(fl.bam, fl.vcf, bs, Homo.sapiens) ## default ideo = FALSE, turned on
## tks <- tracks(fl.bam, fl.vcf, bs, Homo.sapiens, ideo = TRUE)
## tks + xlim(gr17)
gr17 <- GRanges("chr17", IRanges(41234415, 41234569))
tks <- tracks(p.ideo, mismatch = p.mis, dbSNP = p.vr, ref = p.bg, gene = p.txdb,
              heights = c(2, 3, 3, 1, 4)) + xlim(gr17) + theme_tracks_sunset()
png(filename = file.path(outpath, "ggbio.test22.png"), width = 1500, height = 800)
tks
dev.off()


### Buidling circular plot layer by layer
data("CRC", package = "biovizBase")
png(filename = file.path(outpath, "ggbio.test23.png"), width = 1500, height = 800)
autoplot(hg19sub, layout = "circle", fill = "gray70")
dev.off()

p <- ggbio() + circle(hg19sub, geom = "ideo", fill = "gray70") +
  circle(hg19sub, geom = "scale", size = 2) +
  circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
png(filename = file.path(outpath, "ggbio.test24.png"), width = 1500, height = 1500)
p
dev.off()

p <- ggbio() + circle(mut.gr, geom = "rect", color = "steelblue") +
  circle(hg19sub, geom = "ideo", fill = "gray70") +
  circle(hg19sub, geom = "scale", size = 2) +
  circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
png(filename = file.path(outpath, "ggbio.test25.png"), width = 1500, height = 1500)
p
dev.off()


### grandlinear plots
snp <- read.table(system.file("extdata", "plink.assoc.sub.txt", package = "biovizBase"),
                  header = TRUE)
gr.snp <- transformDfToGr(snp, seqnames = "CHR", start = "BP", width = 1)
gr.snp <- keepSeqlevels(gr.snp, as.character(1:22))
data(ideoCyto, package = "biovizBase")
seqlengths(gr.snp) <- as.numeric(seqlengths(ideoCyto$hg18)[1:22])
gr.snp <- gr.snp[!is.na(gr.snp$P)]
values(gr.snp)$pvalue <- -log10(values(gr.snp)$P)

png(filename = file.path(outpath, "ggbio.test26.png"), width = 1500, height = 1500)
autoplot(gr.snp, geom = "point", coord = "genome", aes(y = pvalue))
dev.off()

png(filename = file.path(outpath, "ggbio.test27.png"), width = 1500, height = 1500)
plotGrandLinear(gr.snp, aes(y = pvalue), color = c("#7fc97f", "#fdc086"))
dev.off()

png(filename = file.path(outpath, "ggbio.test27.png"), width = 800, height = 800)
plotGrandLinear(gr.snp, aes(y = pvalue), color = c("#7fc97f", "#fdc086"),
                cutoff = 3, cutoff.color = "blue", cutoff.size = 0.2)
dev.off()

png(filename = file.path(outpath, "ggbio.test28.png"), width = 800, height = 800)
plotGrandLinear(gr.snp, aes(y = pvalue, color = OR), spaceline = TRUE, legend = TRUE)
dev.off()

gro <- GRanges(c("1", "11"), IRanges(c(100, 2e6), width = 5e7))
names(gro) <- c("group1", "group2")
png(filename = file.path(outpath, "ggbio.test29.png"), width = 800, height = 800)
plotGrandLinear(gr.snp, aes(y = pvalue), highlight.gr = gro)
dev.off()


### karyotpye
data(ideoCyto, package = "biovizBase")
png(filename = file.path(outpath, "ggbio.test30.png"), width = 800, height = 800)
autoplot(seqinfo(ideoCyto$hg19), layout = "karyogram")
dev.off()

png(filename = file.path(outpath, "ggbio.test31.png"), width = 800, height = 800)
autoplot(ideoCyto$hg19, layout = "karyogram", cytoband = TRUE)
dev.off()

data(darned_hg19_subset500, package = "biovizBase")
dn <- darned_hg19_subset500
seqlengths(dn) <- seqlengths(ideoCyto$hg19)[names(seqlengths(dn))]
dn <- keepSeqlevels(dn, paste0("chr", c(1:22, "X")))
png(filename = file.path(outpath, "ggbio.test32.png"), width = 800, height = 800)
autoplot(dn, layout = "karyogram")
dev.off()

png(filename = file.path(outpath, "ggbio.test33.png"), width = 800, height = 800)
autoplot(dn, layout = "karyogram", aes(color = exReg, fill = exReg), alpha = 0.5) +
  scale_color_discrete(na.value = "brown")
dev.off()

dn.nona <- dn[!is.na(dn$exReg)]
dn.nona$levels <- as.numeric(factor(dn.nona$exReg))
p.ylim <- autoplot(dn.nona, layout = "karyogram", aes(color = exReg, fill = exReg,
                                                      ymin = (levels - 1) * 10/3,
                                                      ymax = levels * 10 /3))
png(filename = file.path(outpath, "ggbio.test34.png"), width = 800, height = 800)
p.ylim
dev.off()

### Link ranges to your data
data(genesymbol, package = "biovizBase")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
model <- exonsBy(txdb, by = "tx")
model17 <- subsetByOverlaps(model, genesymbol["RBM17"])
exons <- exons(txdb)
exon17 <- subsetByOverlaps(exons, genesymbol["RBM17"])
exon.new <- reduce(exon17)
values(exon.new)$sample1 <- rnorm(length(exon.new), 10, 3)
values(exon.new)$sample2 <- rnorm(length(exon.new), 10, 10)
values(exon.new)$score <- rnorm(length(exon.new))
values(exon.new)$significant <- sample(c(TRUE,FALSE), size = length(exon.new),replace = TRUE)

p17 <- autoplot(txdb, genesymbol["RBM17"])
png(filename = file.path(outpath, "ggbio.test35.png"), width = 800, height = 800)
p17
dev.off()

png(filename = file.path(outpath, "ggbio.test36.png"), width = 800, height = 800)
plotRangesLinkedToData(exon.new, stat.y = c("sample1", "sample2"), annotation = list(p17))
dev.off()


### themes
library(GenomicRanges)
set.seed(1)
N <- 100
gr <- GRanges(seqnames = sample(c("chr1", "chr2", "chr3"),
                                size = N, replace = TRUE),
              IRanges(start = sample(1:300, size = N, replace = TRUE),
                      width = sample(70:75, size = N,replace = TRUE)),
              strand = sample(c("+", "-"), size = N, replace = TRUE),
              value = rnorm(N, 10, 3), score = rnorm(N, 100, 30),
              sample = sample(c("Normal", "Tumor"),
                              size = N, replace = TRUE),
              pair = sample(letters, size = N,
                            replace = TRUE))
seqlengths(gr) <- c(400, 1000, 500)

png(filename = file.path(outpath, "ggbio.test37.png"), width = 800, height = 800)
autoplot(gr)
dev.off()

png(filename = file.path(outpath, "ggbio.test38.png"), width = 800, height = 800)
autoplot(gr) + theme_genome()
dev.off()



