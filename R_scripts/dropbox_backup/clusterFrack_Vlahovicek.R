# gff <- read.delim("final.refined.clusters.no.smr.no.ssi.20151127.proper.gff", header = FALSE, stringsAsFactors = FALSE)
# spl <- do.call(rbind, strsplit(as.character(gff$V9), ";"))
# spl2 <- do.call(rbind, strsplit(spl[,2], "="))[,2]
# 
# spl <- lapply(strsplit(as.character(gff$V9), ";"), function(x) {
#     y <- unlist(strsplit(x, "="))
#     if(length(y) == 2 ) {
#       c(y[2],y[2])
#     } else {
#       c(y[2], y[4])
#     }
# #    length(y)
#   })
# spl <- do.call(rbind, spl)
# 
# gff.aug <- cbind(gff, spl)
# 
# 
# 
# ttt <- 
# by(gff.aug, gff.aug$`2`, function(cl) {
#   cl$`1` <- as.character(cl$`1`)
#   f <- cl[cl$V3 == "exon",]
#   y <- cl[cl$V3 != "exon",]
#   z <- by(f, f$`1`, function(x) {
# #    id <- strsplit(as.character(x[1,10]), "=")[[1]][2]
#     l <- data.frame(x[1,1], x[1,2], "transcript", min(x[,4]), max(x[,5]), x[1,6], x[1,7], x[1,8], 
#                     paste("ID=", x[1,10], ";Parent=", x[1,11], sep=""))
#     yy <- x[,1:9]
#     lab <- paste("ID=", x[1,10], ".", 1:nrow(x), ";Parent=", x[1,10], sep="")
#     yy[,9] <- lab
#     names(l) <- names(yy)
#     rbind(l, yy)
#   })
#   zz <- do.call(rbind, z)
#   rbind(y[,1:9], zz)
#   #  rbind(y, z)
# })
# ttt <- do.call(rbind, ttt)
# head(ttt,20)
# 
# cat("##gff-version 3\n", file="reformatted.gff")
# write.table(ttt, file="reformatted.gff", sep = "\t", append = TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)

library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(Rsamtools)
library(ggbio)
library(biovizBase)
library(pheatmap)

data(ideoCyto)

setwd("C:/Users/kristian/Dropbox/Petr/MedILS2015/Clusters")

calcMask <- function(){
  exogaps <- GRangesList(
    sapply(lncRNAS, function(x) {
      slx <- unique(as.character(seqnames(x)))
      xx <- keepSeqlevels(x, slx)
#      strand(xx) <- "*"
      rx <- range(xx)
      rnas <- subsetByOverlaps(mm9ex, rx, ignore.strand=T)

      # Petr says keep only antisense for the mask. OK.
      rnas <- rnas[strand(rnas) != strand(xx)]
      
      if(length(rnas) > 0) {

  
        # first determine the "coding transcripts 
        # and then filter exons based on whether they belong to a coding transcript
        
        tx <- unique(unlist(rnas$tx_id))
        ctrans <- rep(0, times = length(tx))
        names(ctrans) <- as.character(tx)
        for(i in 1:length(rnas)) {
          cds <- unlist(rnas$cds_id[i])
          txe <- unlist(rnas$tx_id[i])
          if(length(cds) == length(txe))
            ctrans[as.character(txe)] <- ctrans[as.character(txe)] + !is.na(cds)
        }
        ctrEx <- sapply(rnas$tx_id, function(x) sum(ctrans[as.character(unlist(x))])) > 0
        simplEx <- sapply(rnas$cds_id, function(x) sum(!is.na(x)) > 0 )
        
        codingEx <- simplEx | ctrEx
        
        fovl <- findOverlaps(rnas[codingEx == FALSE], rnas[codingEx])
        if(length(fovl) > 0) {
          codingEx[which(codingEx == FALSE)[unique(queryHits(fovl))]] <- TRUE
        }
        rnas <- rnas[codingEx]
      }
      #    if(length(rnas) > 0) {
      #      codingEx <- tapply(rnas$cds_id, unlist(rnas$gene_id), function(x) {
      #        ulid <- unlist(x)
      #        coding <- sum(!is.na(ulid)) > 0
      #        rep(coding, times = length(x))
      #      })
      #      codingEx <- unlist(codingEx)
      #     }
      
      
      if(length(rnas) > 0) {
        rnas <- keepSeqlevels(rnas, unique(c(as.character(seqnames(rnas)), slx)))
        #      cat(unique(c(as.character(seqnames(rnas)), slx)), "\n")
        strand(rnas) <- "*"
        rl <- range(rnas)
        #      cat(width(rnas), "\n")
        strand(rl) <- "*"
        tx <- range(c(rl, rx), ignore.strand = T)
        gx <- gaps(rnas, start = start(tx), end = end(tx))
        #    gx <- trim(gaps(xx))
        gx <- gx[strand(gx) == "*"]
        subsetByOverlaps(gx, rx)
      } else {
        rx
      }
    })
  )
  zzz <- sapply(exogaps, function(x) length(x) < 1)
  exogaps <- exogaps[!zzz]
  exogaps
}

mastertable <- read.delim("../ngs.est.ss.rrna.repeats.junction.fpkm.strict.one.20160204.retro.cpat.muerv.novelty.final.txt",
                          stringsAsFactors = FALSE)
mastertable$K_ID <- paste("cluster_", sub("gv_2c_bl_", "", mastertable$cluster_ID), sep="")
mastertable$simple_stage <- factor(c("ZY", "ZY", "EM", "MA", "MA", "EM")[factor(mastertable$max_FPKM_stage)], 
                                   levels=c("MA", "ZY", "EM"))
mastertable$simple_stage_levels <- as.numeric(mastertable$simple_stage)
rownames(mastertable) <- mastertable$K_ID 

master.gr <- makeGRangesFromDataFrame(df = mastertable, keep.extra.columns = TRUE)
seqlengths(master.gr)<-seqlengths(ideoCyto$mm9)[names(seqlengths(master.gr))]

mm9UCSC <- makeTxDbFromUCSC(genome="mm9", tablename="knownGene")
mm9ENS  <- makeTxDbFromUCSC(genome="mm9", tablename="ensGene")
cols <- c("gene_id", "tx_id", "exon_id", "cds_id")
mm9ex <- c(exons(mm9UCSC, columns=cols), exons(mm9ENS, columns=cols))

chrfile <- file.path(".","mm9.chrlen.txt")
chrominfo <- read.delim(chrfile, header = FALSE, stringsAsFactors = FALSE)
names(chrominfo) <- c("chrom", "length")
chrominfo$is_circular <- FALSE
chrominfo[chrominfo$chrom == "chrM", "is_circular"] <- TRUE

sInfo <- Seqinfo(seqnames   = chrominfo$chrom, 
                 seqlengths = chrominfo$length,
                 isCircular = chrominfo$is_circular,
                 genome     = "mm9")

# (txdb <- makeTxDbFromGFF("reformatted.gff", format="gff", 
#                          chrominfo = chrominfo,
#                          organism = "Mus musculus"))
(txdb <- makeTxDbFromGFF("http://hex.bioinfo.hr/~rosa/FPKM_to_check/final.refined.clusters.no.smr.no.ssi.20160316.1600.gff", format="gff", 
#(txdb <- makeTxDbFromGFF("http://hex.bioinfo.hr/~rosa/FPKM_to_check/final.refined.clusters.no.smr.no.ssi.20160307.1602.gff", format="gff", 
#(txdb <- makeTxDbFromGFF("final.refined.clusters.no.smr.no.ssi.20160307.1602.gff", format="gff", 
                         chrominfo = chrominfo,
                         organism = "Mus musculus"))

lncRNAS <- genes(txdb)
unlistedReducedExons <- unlist(reduce(exonsBy(txdb, "gene")))
unlistedReducedExons$gene_id <- names(unlistedReducedExons)
lncTranscripts <- transcripts(txdb, columns = c("gene_id", "tx_name"))

exogaps <- calcMask()

maskReads <- function(reads, mask = exogaps, limits = lncRNAS) {
  r <- subsetByOverlaps(reads, limits, type="within", ignore.strand = T)
  subsetByOverlaps(r, mask, ignore.strand = T)
}

countJunctions <- function(features, reads, ignore.strand=T, inter.feature=T) {
  fovl <- findOverlaps(features, reads, ignore.strand=ignore.strand)
  fosp <- split(njunc(reads[subjectHits(fovl)]), queryHits(fovl))
  sJunc <- sapply(fosp, sum)
  res <- setNames(rep(0, times = length(features)), names(features))
  #  names(res) <- names(features)
  res[as.integer(names(sJunc))] <- sJunc
  res
}

countJunctionsInBam <- function(f, features, param = p1, mask = exogaps, limits = lncRNAS) {
  reads<- readGappedReads(f, param = param)
  reads.f <- maskReads(reads, mask = mask, limits = limits)
  
  countJunctions(features, reads.f)
  
}

p1 <- ScanBamParam(which=unlist(exogaps))

(ebg <- exonsBy(txdb, by="gene"))
seqinfo(ebg) <-  sInfo
djWidths <- sapply(disjoin(ebg), function(x) sum(width(x)))

(ebt <- exonsBy(txdb, by="tx", use.names=TRUE))

infiles <- system2("find", "/common/DB/vfranke/Base/Encode/mm9/RNASeq/Mapped/ -name \\*.bam", stdout=TRUE)
bamfiles <- BamFileList(infiles, yieldSize=200000, asMates=FALSE)
bamfiles1 <- BamFileList(infiles[1], yieldSize=200000, asMates=FALSE)

#se <- summarizeOverlaps(features=ebg, reads=bamfiles,
#                        mode="IntersectionNotEmpty",
se <- summarizeOverlaps(features=unlistedReducedExons, reads=bamfiles,
                        mode="Union",
                        inter.feature=FALSE,
                        singleEnd=TRUE,
                        ignore.strand=TRUE, 
                        fragments = FALSE,
                        param=p1, 
                        preprocess.reads=maskReads)

#save(se, file="count_se.Robj")
save(se, file="count_reduced_masked_encode_if.Robj")

se.junc <- sapply(bamfiles, function(f) {
  countJunctionsInBam(f, ebg)
})

#save(se, file="count_se.Robj")
save(se.junc, file="count_junctions_exonsbygenes_masked_encode_if.Robj")

#ex <- exons(txdb, columns=c("gene_id", "tx_name", "exon_name"))
#seqinfo(ex) <-  sInfo[seqlevels(ex)]


#se.ex <- summarizeOverlaps(features=ex, reads=bamfiles,
#                        mode="Union",
#                        singleEnd=FALSE,
#                        ignore.strand=TRUE)

#totalMdata <- cbind(mcols(ex), as.data.frame(assay(se.ex)))
#names(totalMdata) <- sub(".bam", "", names(totalMdata))
#mcols(ex) <- totalMdata

#save(se.ex, file="count_se_ex.Robj")
#save(ex, file="exon_counts.Robj")


infiles.fu <- system2("find", "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd/Star_PairEnd_EnsemblAnnot -name \\*.bam | grep -v uniq", stdout=TRUE)
bamfiles.fu <- BamFileList(infiles.fu, yieldSize=200000, asMates=FALSE)

#se.ex.fu <- summarizeOverlaps(features=ex, reads=bamfiles.fu,
se.ex.fu <- summarizeOverlaps(features=unlistedReducedExons, reads=bamfiles.fu,
                              mode="Union",
                              inter.feature=FALSE,
                              singleEnd=TRUE,
                              ignore.strand=TRUE,
                              fragments = FALSE,
                              param=p1, 
                              preprocess.reads=maskReads)

#totalMdata.fu <- cbind(mcols(ex), as.data.frame(assay(se.ex.fu)))
#names(totalMdata.fu) <- sub(".bam", "", names(totalMdata.fu))
#mcols(ex) <- totalMdata.fu


#save(se.ex.fu, file="count_se_ex_fu.Robj")
save(se.ex.fu, file="count_reduced_masked_fu_if.Robj")
#save(ex, file="exon_counts_fugaku_tissue.Robj")

se.fu.junc <- sapply(bamfiles.fu, function(f) {
  countJunctionsInBam(f, ebg)
})

#save(se, file="count_se.Robj")
save(se.fu.junc, file="count_junctions_exonsbygenes_masked_fu_if.Robj")



infiles.su <- system2("find", "/common/WORK/kristian/Projekti/Petr/LNCRNA/Susor_data/Mapping/bbmap/mm9 -name \\*.bam", stdout=TRUE)
bamfiles.su <- BamFileList(infiles.su, yieldSize=200000, asMates=FALSE)

#se.ex.fu <- summarizeOverlaps(features=ex, reads=bamfiles.fu,
se.ex.su <- summarizeOverlaps(features=unlistedReducedExons, reads=bamfiles.su,
                              mode="Union",
                              inter.feature=FALSE,
                              singleEnd=TRUE,
                              ignore.strand=TRUE,
                              fragments = FALSE,
                              param=p1, 
                              preprocess.reads=maskReads)

#totalMdata.fu <- cbind(mcols(ex), as.data.frame(assay(se.ex.fu)))
#names(totalMdata.fu) <- sub(".bam", "", names(totalMdata.fu))
#mcols(ex) <- totalMdata.fu


#save(se.ex.fu, file="count_se_ex_fu.Robj")
save(se.ex.su, file="count_reduced_masked_susor_if.Robj")

se.su.junc <- sapply(bamfiles.su, function(f) {
  countJunctionsInBam(f, ebg)
})

#save(se, file="count_se.Robj")
save(se.su.junc, file="count_junctions_exonsbygenes_masked_su_if.Robj")

#junctions <- cbind(se.junc, se.fu.junc, se.su.junc)
junctions <- cbind(se.junc, se.fu.junc)
colnames(junctions) <- paste(colnames(junctions), "junc", sep=".")

save(junctions, file="count_junctions_bygenes_masked_encode_fugaku_susor_if.Robj")


totalMdata <- cbind(mcols(unlistedReducedExons), 
                    as.data.frame(assay(se)),
                    as.data.frame(assay(se.ex.fu)),
                    as.data.frame(assay(se.ex.su)))
                    

names(totalMdata) <- sub(".bam", "", names(totalMdata))
uRexFinal <- unlistedReducedExons
mcols(uRexFinal) <- totalMdata

save(uRexFinal, file="unlisted_reduced_exons_lncRNA_counts_encode_fugaku_susor_if_se.Robj")


## Windows

#load("count_se.Robj")
#rawCounts <- assay(se)
#nrawCounts <- rawCounts / djWidths[rownames(rawCounts)]

tissueList <- c("CnsE11half", "CnsE14", "CnsE18", "WbrainE14half", "FlobeAdult8wks", "CortexAdult8wks",
                "CbellumAdult8wks", "StomAdult8wks", "LiverE14", "LiverE18", "LiverAdult8wks",
                "DuodAdult8wks", "SmintAdult8wks", "LgintAdult8wks", "ColonAdult8wks", "LimbE14half",
                "LungAdult8wks", "HeartAdult8wks", "BladderAdult8wks", "KidneyAdult8wks", "ThymusAdult8wks",
                "MamgAdult8wks", "SpleenAdult8wks", "SfatAdult8wks", "GfatAdult8wks", "OvaryAdult8wks",
                "TestisAdult8wks","PlacAdult8wks")

tissueList <- c("CnsE11half", "CnsE14", "CnsE18", "FlobeAdult8wks", "CortexAdult8wks",
                "CbellumAdult8wks", "StomAdult8wks", "LiverAdult8wks",
                "DuodAdult8wks", "SmintAdult8wks", "LgintAdult8wks", "ColonAdult8wks", 
                "LungAdult8wks", "HeartAdult8wks", "BladderAdult8wks", "KidneyAdult8wks", "ThymusAdult8wks",
                "MamgAdult8wks", "SpleenAdult8wks", "OvaryAdult8wks",
                "TestisAdult8wks","PlacAdult8wks")

# GV-MII-1C-1CDNAm-2C-2CDNAm-4C-Mo-Bl-MIIPA-1CPA
tissueList <- c("s_GV.WE", "s_MII.WE", "s_1cell.WE", "s_1cell.WE_DNAm", 
                "s_2cell.WE", "s_2cell.WE_DNAm", "s_4cell.WE", 
                "s_Morula.WE", "s_Blast.WE", "s_MII.PA", "s_1cell.PA")

# no inhibitors
tissueList <- c("s_GV.WE", "s_MII.WE", "s_1cell.WE", 
                "s_2cell.WE", "s_4cell.WE", 
                "s_Morula.WE", "s_Blast.WE", "s_MII.PA", "s_1cell.PA")

# no inhibitors, no PA
tissueList <- c("s_GV.WE", "s_MII.WE", "s_1cell.WE", 
                "s_2cell.WE", "s_4cell.WE", 
                "s_Morula.WE", "s_Blast.WE")


# GV-2C-BL
tissueList <- c("s_GV.WE", "s_2cell.WE", "s_Blast.WE")

# Susor data

tissueList <- c("X7510_sorted", "X7511_sorted", "X7512_sorted", "X7877_sorted")

# All (relevant)

tissueList <- c("CnsE11half", "CnsE14", "CnsE18", "FlobeAdult8wks", "CortexAdult8wks",
                "CbellumAdult8wks", "StomAdult8wks", "LiverAdult8wks",
                "DuodAdult8wks", "SmintAdult8wks", "LgintAdult8wks", "ColonAdult8wks", 
                "LungAdult8wks", "HeartAdult8wks", "BladderAdult8wks", "KidneyAdult8wks", "ThymusAdult8wks",
                "MamgAdult8wks", "SpleenAdult8wks", "OvaryAdult8wks",
                "TestisAdult8wks","PlacAdult8wks",
                "s_GV.WE", "s_MII.WE", "s_1cell.WE", 
                "s_2cell.WE", "s_4cell.WE", 
                "s_Morula.WE", "s_Blast.WE", "s_MII.PA", "s_1cell.PA",
                "X7510_sorted", "X7511_sorted", "X7512_sorted", "X7877_sorted")



load("exon_counts_fugaku_tissue.Robj")
load("unlisted_reduced_exons_lncRNA_counts_encode_fugaku_susor.Robj")
load("unlisted_reduced_exons_lncRNA_counts_encode_fugaku_susor_if_se.Robj")


# remove antisense overlapping exons for expression analysis!
tex <- ex[unlist(ex$gene_id) %in% mastertable$K_ID]

kgOver <- findOverlaps(mm9ex, tex, ignore.strand=TRUE)
sstrand <- strand(tex[subjectHits(kgOver)])
qstrand <- strand(mm9ex[queryHits(kgOver)])
exclude <- subjectHits(kgOver[sstrand != qstrand])

cex <- tex[-exclude]

cex <- uRexFinal

libsizefile <- file.path(".","libsizes.txt")
libSizes <- read.delim(libsizefile, header = FALSE, stringsAsFactors = FALSE)
names(libSizes) <- c("track", "size")
rownames(libSizes) <- libSizes$track

libSizes$size <- libSizes$size / 1e6
geneWidth <- width(cex) / 1000
names(geneWidth) <- cex$exon_name
names(geneWidth) <- make.names(cex$gene_id, unique = TRUE)

mData <- mcols(cex)
rownames(mData) <- mData$exon_name
rownames(mData) <- make.names(mData$gene_id, unique = TRUE)
validTracks <- base::intersect(tissueList, base::intersect(names(mData), libSizes$track))

libnorm <- sweep(as.matrix(mData[,validTracks]),2,libSizes[validTracks, "size"], '/')
RPKM    <- sweep(libnorm, 1, geneWidth[row.names(libnorm)], '/')

save(RPKM, file="allRPKM.Robj")

maxRPKM <- aggregate(as.data.frame(RPKM), list(eid = unlist(mData$gene_id)), FUN = max)
rownames(maxRPKM) <- maxRPKM$eid
maxRPKM$eid <- NULL

bigguys <- apply(maxRPKM, 1, function(x) any(x > 1))
bigguys4 <- apply(maxRPKM, 1, function(x) any(x > 4))
#bigguys <- apply(maxRPKM, 1, function(x) any(x > 0))
# for all clusters:
# bigguys <- rep(TRUE, nrow(maxRPKM))
  
annoCol <- data.frame(
#  chromosome = sub("cluster_(chr\\w+)_\\d+","\\1",rownames(maxRPKM)),
                      dev_stage = mastertable[rownames(maxRPKM), "simple_stage"])
rownames(annoCol) <- rownames(maxRPKM)

maxInClust <- apply(maxRPKM, 1, max)
maxRPKMsc <- apply(maxRPKM, 2, function(x) x/maxInClust)

pheatmap(t(log(maxRPKM[bigguys,]+2)), cluster_rows = FALSE, annotation_col = annoCol)
pheatmap(t(log(maxRPKM[bigguys,]+2)), cluster_rows = FALSE, show_colnames = FALSE)

pheatmap(t(log(maxRPKMsc[bigguys,]+2)), cluster_rows = FALSE, cutree_cols = 6)
zz <- pheatmap(t(log(maxRPKMsc[bigguys,]+2)), cluster_rows = FALSE, cutree_cols = 6,color = colorRampPalette(c("white", "black"))(100))

pheatmap(t(log(maxRPKM[bigguys,]+2)), cluster_rows = FALSE, color = colorRampPalette(c("black", "yellow" , "orange", "red"))(100))
pheatmap(t(log(maxRPKM[bigguys,]+2)), cluster_rows = FALSE, color = colorRampPalette(c("white", "black"))(100))


classes <- cutree(zz$tree_col, 6)
# constraints

mRm <- as.data.frame(maxRPKMsc[bigguys,])

# for no inhibitor GV-Bl, slide 22 Rehovot:

mRmNoMo <- as.data.frame(maxRPKM[bigguys,])
mRmNoMo <- mRmNoMo[mRmNoMo$s_Morula.WE == 0 & mRmNoMo$s_GV.WE != 0,]
pheatmap(t(log(mRmNoMo[order(mRmNoMo$s_GV.WE/mRmNoMo$s_2cell.WE),]+2)), 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         color = colorRampPalette(c("black", "yellow" , "orange", "red"))(100))

pheatmap(t(log(mRmNoMo[order(mRmNoMo$s_MII.WE, 
                             mRmNoMo$s_GV.WE, 
                             
                             mRmNoMo$s_1cell.WE, 
                             mRmNoMo$s_2cell.WE, 
                             mRmNoMo$s_4cell.WE, 
                             mRmNoMo$s_Morula.WE, 
                             mRmNoMo$s_Blast.WE),]+2)), 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         color = colorRampPalette(c("black", "yellow" , "orange", "red"))(100))

classes <- vector()
classes[row.names(mRm)] <- NA

#manual stratification
classes[mRm$s_Blast.WE == 1 & mRm$s_GV.WE <= mRm$s_2cell.WE] <- 1
classes[mRm$s_Blast.WE == 1 & mRm$s_GV.WE >  mRm$s_2cell.WE] <- 2

classes[mRm$s_2cell.WE == 1 & mRm$s_GV.WE <= mRm$s_Blast.WE] <- 3
classes[mRm$s_2cell.WE == 1 & mRm$s_GV.WE >  mRm$s_Blast.WE] <- 4

classes[mRm$s_GV.WE == 1 & mRm$s_2cell.WE <= mRm$s_Blast.WE] <- 5
classes[mRm$s_GV.WE == 1 & mRm$s_2cell.WE >  mRm$s_Blast.WE] <- 6
# re-class the ones with 2C = BL = 0 to 6
classes[mRm$s_GV.WE == 1 & mRm$s_2cell.WE == 0 & mRm$s_Blast.WE == 0] <- 6

mRmOrd <- mRm[order(classes, mRm$s_GV.WE, mRm$s_2cell.WE, mRm$s_Blast.WE),]
gaps <- cumsum(rle(classes[order(classes)])$lengths)

pheatmap(t(log(mRmOrd+2)), 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         gaps_col = gaps,
         color = colorRampPalette(c("white", "black"))(100))


maxRPKMmelt <- melt(as.matrix(mRm), varnames = c("cluster", "stage"))
maxRPKMmelt$stage <- factor(maxRPKMmelt$stage, levels = c("s_GV.WE", "s_2cell.WE", "s_Blast.WE"), 
                            labels = c("GV", "2C", "BL"))
maxRPKMmelt$panel <- factor(classes[maxRPKMmelt$cluster], levels = c(6,4,2,5,3,1))

ggplot(maxRPKMmelt, aes(x=stage, y=value, group=cluster)) + 
  geom_line(colour="cornflowerblue", size=1) + 
#  geom_smooth(aes(group=panel), method="loess", se=FALSE, colour="red", size=2) +
  stat_summary(aes(group=panel), fun.y = mean, fun.ymin = min, fun.ymax = max,
               colour = "red", geom="line", size=2) +
  facet_wrap(~panel) +
  coord_fixed(ratio=2, ylim = c(0, 1), xlim = c(1, 3)) +
  scale_y_continuous(breaks = c(0, .5, 1), limits = c(0, 1)) +
  theme_bw() +
  theme(panel.border = element_rect(linetype = "solid", colour = "black", size = 1),
        strip.background = element_blank(), strip.text = element_blank(),
        panel.margin.x = unit(0.7, "lines"),
        panel.grid.major = element_line(size=1), 
        axis.line = element_blank(),
        axis.ticks = element_blank())


ggplot(maxRPKMmelt, aes(x=stage, y=value, group=cluster)) + 
  geom_line(colour="black", size=.8, alpha = .3) + 
  #  geom_smooth(aes(group=panel), method="loess", se=FALSE, colour="red", size=2) +
#  stat_summary(aes(group=panel), fun.y = mean, fun.ymin = min, fun.ymax = max,
#               colour = "red", geom="line", size=2) +
#  facet_wrap(~panel) +
  coord_fixed(ratio=2, ylim = c(0, 1), xlim = c(1, 3)) +
  scale_y_continuous(breaks = c(0, .5, 1), limits = c(0, 1)) +
  theme_bw() +
  theme(panel.border = element_rect(linetype = "solid", colour = "black", size = 1),
        strip.background = element_blank(), strip.text = element_blank(),
        panel.margin.x = unit(0.7, "lines"),
        panel.grid.major = element_line(size=1), 
        axis.line = element_blank(),
        axis.ticks = element_blank())

# MT clusters
MTs <- list(MTA = c("MTA"),
            MTB = c("MTB"),
            MTC = c("MTC"), 
            MTD = c("MTD"), 
            MTE = c("MTE"),
            MT2 = c("MT2_Mm"),
            MTAtoD = c("MTA|MTB|MTC|MTD"))

mtClusters <- sapply(names(MTs), function(x){
  rn <- row.names(mastertable[grep(MTs[x], mastertable$over_TSS_and_promoter_50bp),])
  rnr <- intersect(rn, row.names(maxRPKM))
  matr <- maxRPKM[rnr, ]+2
  if(x == "MT2_Mm") {
    matr <- matr[order(matr$s_2cell.WE), ]
  } else {
    matr <- matr[order(matr$s_GV.WE), ]
  }
  
  pheatmap(t(log(matr)), 
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           color = colorRampPalette(c("black", "yellow" , "orange", "red"))(100),
           main = x,
           filename = paste(x, "pdf", sep="."),
           width = 15, 
           height = 10)
  rnr
})



# plot karyogram

autoplot(master.gr, layout="karyogram", aes(color=simple_stage, fill =simple_stage)) + 
  scale_colour_brewer(palette="Set1") + 
  scale_fill_brewer(palette="Set1") +
  theme_clear()

autoplot(master.gr, layout="karyogram", aes(color=simple_stage, fill =simple_stage, ymin = (simple_stage_levels - 1) *10/3, ymax = simple_stage_levels *10 /3 ))


# Range fddling for quantification

lncRNAS <- genes(txdb)
hits <- findOverlaps(lncRNAS, mm9ex, ignore.strand=TRUE)
mm9subset <- mm9ex[subjectHits()]

sapply(lncRNAS, function(x) {
  fo <- findOverlaps(x, mm9subset)
  ovlEx <- mm9subset[subjectHits(fo), drop = TRUE]
  neg <- reduce(gaps(ovlEx, start = start(x), end = end(x)), drop=TRUE)
  neg
})


mm9subset <- mm9ex[subjectHits(findOverlaps(lncRNAS,mm9ex, ignore.strand=TRUE))]
p1 <- ScanBamParam(which=lncRNAS)
gapped_alignments <- readGappedReads("chr4.part.bam", use.names = T, param = p1)
gapped_subset <- subsetByOverlaps(gapped_alignments, gaps(mm9subset))


tt <- lncRNAS[lncRNAS$gene_id == "cluster_chr4_594"]
tt <- lncRNAS[lncRNAS$gene_id == "cluster_chr4_335"]
tt <- lncRNAS[lncRNAS$gene_id == "cluster_chr4_380"]
tt <- lncRNAS[lncRNAS$gene_id == "cluster_chr1_11"]
tt <- lncRNAS[lncRNAS$gene_id == "cluster_chr10_397"]
tt <- lncRNAS[lncRNAS$gene_id == "cluster_chr8_76"]
tt <- lncRNAS[lncRNAS$gene_id == "cluster_chr8_337"]
tt <- lncRNAS[lncRNAS$gene_id == "cluster_chr8_336"]
tt <- lncRNAS[lncRNAS$gene_id == "cluster_chrX_482"]
tt <- lncRNAS[lncRNAS$gene_id == "cluster_chrX_445"]

mrnas <- subsetByOverlaps(mm9ex, tt, ignore.strand=TRUE)
lnc <- reduce(ex[unlist(ex$gene_id) == "cluster_chr4_594"])
p1 <- ScanBamParam(which=mrnas[3:7])
#gapped_alignments <- readGappedReads("chr4.part.bam", use.names = T, param = p1)
gapped_alignments <- readGappedReads("lncRNA_filtered_GV.bam", use.names = T, param = p1)
gapped_alignments <- readGappedReads("lncRNA_filtered_Testis.bam", use.names = T, param = p1)

grt <- GeneRegionTrack(txdb, stacking="dense", geneSymbols=TRUE, name = "lncRNAs", showFeatureId=TRUE)
grt <- GeneRegionTrack(txdb, stacking="squish", geneSymbols=TRUE, name = "lncRNAs", showFeatureId=TRUE)

allIntrons <- intronsByTranscript(mm9UCSC)
myIntrons <- subsetByOverlaps(allIntrons, tt, ignore.strand=T)
grtInt <- AnnotationTrack(disjoin(unlist(myIntrons)), name = "Introns")

grtUC <- GeneRegionTrack(mm9UCSC, fill="blue", name = "UCSC", showFeatureId=TRUE)
grtEN <- GeneRegionTrack(mm9ENS, fill="red", name = "Ensembl", showFeatureId=TRUE)
#alTrack <- AlignmentsTrack(range="chr4.part.bam", isPaired=FALSE)
alTrack <- AlignmentsTrack(range="lncRNA_filtered_GV.bam", 
                           isPaired=FALSE, name="AllReads", col.gap = "orange", col.mates = "purple")
alTrack <- AlignmentsTrack(range="lncRNA_filtered_Testis.bam", 
                           isPaired=FALSE, name="AllReads", col.gap = "orange", col.mates = "purple")

exonicus <- subsetByOverlaps(exonsBy(mm9UCSC, by="gene"), tt, ignore.strand=T)
exonicus <- subsetByOverlaps(mm9ex, lncRNAS, ignore.strand=T)
exogaps <- calcMask()

save(exogaps, file="exogaps.RData")

grtGa <- AnnotationTrack(exogaps, name = "MASK")

pTrack.eg <- function(r){
  rr <- exogaps[[r]]
  rl <- lncRNAS[r]
#  sol <- subsetByOverlaps(gapped_alignments, rr, ignore.strand=T, type="within")
#  sol <- subsetByOverlaps(gapped_alignments, rr, ignore.strand=T)
#  sol <- subsetByOverlaps(sol, rl, ignore.strand=T, type="within")
  
  sol <- maskReads(gapped_alignments, rr, rl)
  
  alT <- as(sol, "GRanges")
  alT$cigar <- cigar(sol)
  alTrackF <- AlignmentsTrack(alT, name="MaskedOut reads", col.gap = "blue", type="pileup")
  
  solEx <- subsetByOverlaps(sol, ebg[[r]], ignore.strand=T)
  solExTr <- as(solEx,  "GRanges")
  solExTr$cigar <- cigar(solEx)
  exTracks <- AlignmentsTrack(solExTr, name="Final reads", col.gap = "green", type="pileup")

  juncTracks <- GeneRegionTrack(summarizeJunctions(solEx), name="Junctions", fill="lightgreen")
  
  st <- start(range(rr))
  en <- end(range(rr))
  ch <- as.character(seqnames(rr))
  id <- r
  ti <- paste(id, " ", ch, ":", st, "-", en, sep = "")
  
  
  plotTracks(c(alTrack, alTrackF, exTracks, juncTracks, grtGa, grtUC, grtEN, grt), 
             from = st, to = en, chromosome = ch,
             main = ti)
  
}


plotTracks(c(alTrack, grtInt, grtEx, grtGa, grtUC, grt), from = 140300272, to = 140306640, col.mates = "purple", col.gap = "orange")
plotTracks(c(alTrack, alTrackF, grtGa, grtUC, grtEN, grt), from = 140300272, to = 140306640, chromosome = "chr4", col.mates = "purple", col.gap = "orange")
plotTracks(c(alTrack, alTrackF, grtGa, grtUC, grtEN, grt), from = 131293342, to = 131301710, chromosome = "chr4", col.mates = "purple", col.gap = "orange")
plotTracks(c(alTrack, alTrackF, grtGa, grtUC, grtEN, grt), from = 142086218, to = 142120411, chromosome = "chr4", col.mates = "purple", col.gap = "orange")
plotTracks(c(alTrack, grtGa, grtUC, grtEN, grt), from = 9931412, to = 9934333, chromosome = "chr1", col.mates = "purple", col.gap = "orange")
plotTracks(c(alTrack, grtGa, grtUC, grtEN, grt), from = 64726440, to = 64737219, chromosome = "chr1", col.mates = "purple", col.gap = "orange")
plotTracks(c(alTrack, grtGa, grtUC, grtEN, grt), from = 76485667, to = 76498569, chromosome = "chr1", col.mates = "purple", col.gap = "orange")
plotTracks(c(alTrack, grtGa, grtUC, grtEN, grt), from = 173047305, to = 173055302, chromosome = "chr1", col.mates = "purple", col.gap = "orange")
plotTracks(c(alTrack, grtGa, grtUC, grtEN, grt), from = 110509036, to = 110530061, chromosome = "chr10", col.mates = "purple", col.gap = "orange", geneSymbols=TRUE)
plotTracks(c(alTrack, grtGa, grtUC, grtEN, grt), from = 20019322, to = 20020392, chromosome = "chr8", col.mates = "purple", col.gap = "orange", geneSymbols=TRUE)
plotTracks(c(alTrack, grtGa, grtUC, grtEN, grt), from = 32170348, to = 32174421, chromosome = "chr15", col.mates = "purple", col.gap = "orange", geneSymbols=TRUE)
plotTracks(c(alTrack, grtGa, grtUC, grtEN, grt), from = 23746674, to = 23748757, chromosome = "chrX", col.mates = "purple", col.gap = "orange", geneSymbols=TRUE)
plotTracks(c(alTrack, grtGa, grtUC, grtEN, grt), from = 163500562, to = 163508321, chromosome = "chrX", col.mates = "purple", col.gap = "orange", geneSymbols=TRUE)
plotTracks(c(alTrack, grtGa, grtUC, grtEN, grt), from = 32717279, to = 33164416, chromosome = "chrX", col.mates = "purple", col.gap = "orange", geneSymbols=TRUE)

pdf("masks_uber_max_final_testis.pdf", width = 10, height = 20)
for(i in 1:length(lncRNAS)) {
  st <- start(lncRNAS[i])
  en <- end(lncRNAS[i])
  ch <- as.character(seqnames(lncRNAS[i]))
  id <- lncRNAS[i]$gene_id
  ti <- paste(id, " ", ch, ":", st, "-", en, sep = "")
  
  mask <- exogaps[[id]]
  
#  sol <- subsetByOverlaps(gapped_alignments, mask, ignore.strand=T)
#  sol <- subsetByOverlaps(sol, lncRNAS[i], ignore.strand=T, type="within")
  
  sol <- maskReads(gapped_alignments, mask, lncRNAS[i])
  
  
  alT <- as(sol, "GRanges")
  alT$cigar <- cigar(sol)
  alTrackF <- AlignmentsTrack(alT, name="MaskedOut reads", col.gap = "blue", type="pileup")
  
  solEx <- subsetByOverlaps(sol, ebg[[id]], ignore.strand=T)
  solExTr <- as(solEx,  "GRanges")
  solExTr$cigar <- cigar(solEx)
  exTracks <- AlignmentsTrack(solExTr, name="Final reads", col.gap = "green", type="pileup")

  juncTracks <- GeneRegionTrack(summarizeJunctions(solEx), name="Junctions", col="lightgreen", fill="lightgreen")
  
  
  cat(ti, "\n")
  
  plotTracks(c(alTrack, alTrackF, exTracks, juncTracks, grtGa, grtUC, grtEN, grt), 
             from = st, 
             to = en, 
             chromosome = ch,
             main = ti)
}
dev.off()


install.packages(c('digest', 'Hmisc', 'matrixStats', 'Rcpp', 'RCurl', 'scales', 'XML'))
biocLite(c('graph', 'RBGL', 'rtracklayer', 'VariantAnnotation', 'zlibbioc'))
