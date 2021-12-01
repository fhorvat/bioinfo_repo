library("GenomicRanges")
library("Biostrings")
library("BSgenome.Mmusculus.UCSC.mm10")
library("seqinr")
library("dplyr")
library("ggplot2")

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/")
chr='chr14'
elements_all <- read.delim("MT_MT2_ORR_MLT_allElements.txt", header = T, stringsAsFactors = F)
elements <- elements_all

MT <- elements[grepl("MT", elements$repName), ]
MT <- MT[!grepl("MT2", MT$repName), ]
MT <- MT[!grepl("-int", MT$repName), ]
MT <- makeGRangesFromDataFrame(MT, keep.extra.columns = TRUE)
MT <- MT[seqnames(MT) == chr, ]
MT_start <- data.frame(start = start(MT), class = "MT")

MT2 <- elements[grepl("MT2", elements$repName), ]
MT2 <- MT2[!grepl("-int", MT2$repName), ]
MT2 <- makeGRangesFromDataFrame(MT2, keep.extra.columns = TRUE)
MT2 <- MT2[seqnames(MT2) == chr, ]
MT2_start <- data.frame(start = start(MT2), class = "MT2")

all_start <- full_join(MT_start, MT2_start)

library(BSgenome.Mmusculus.UCSC.mm10)

# ggplot(data = all_start, aes(start, fill = class)) + geom_tile()
# ggplot()+  geom_density(data = all_start, aes(start, color = class), adjust=1/20, size=1.5) 

# library(GEOquery)
# geo=getGEO('GSE56879')
# d=as.character(pData(geo[[1]])$supplementary_file_1)
# d = d[1:26]
# lapply(1:26, function(x)download.file(d[x], destfile=file.path('~/Tmp',basename(d[x]))))
files = list.files('~/Tmp', pattern='txt', full.names=T)
files = files[grepl(chr, files)]
files = files[!grepl('bulk', files)]

library(data.table)
dat = lapply(files, function(x)fread(x))
set = ifelse(grepl('MII', files),'MII','mES')
uset = rbindlist(dat)
uset[,set := rep(set, times=sapply(dat, nrow))]
uset[,exp := rep(1:length(dat), times=sapply(dat, nrow))]

library(AnnotationHub)
library(rtracklayer)
h = AnnotationHub()
chain = h[['AH14596']]
mm9.cpg = h[["AH6117"]]
mm10.cpg = unlist(liftOver(mm9.cpg , chain = chain))
cpg.ind = countOverlaps(with(uset, GRanges(paste0('chr',V1), IRanges(V2,width=1))), mm9.cpg) >0
uset = uset[!cpg.ind,]


br = seq(1,seqlengths(Mmusculus)[chr],by=1e6)
bins = cut(all_start$start, 
           breaks=br,
           labels= as.character(head(br,-1)))
all_start = all_start[order(as.numeric(all_start$start)),]
all_start$bin=bins
all_start = all_start[order(as.numeric(all_start$start)),]

bins2 = cut(uset$V2, 
            breaks=br,
            labels= as.character(head(br,-1)))
uset[, bin := bins2]
uset$perc = with(uset, V4/(V4+V5+1))
# uset[, perc := V4/(V4+V5+1)]
# metcov = uset[,list(perc=mean(perc), cov=sum(V4)), by=list(bin, set)]
metsum = dcast(data=uset, formula=bin~set, value.var='V4', fill=0, fun.aggregate='sum')

metcov = dcast(data=uset, formula=bin~set+exp, value.var='perc', fill=0, fun.aggregate='mean')
metcov = metcov[with(metsum, MII>100 | mES > 100),]
metcov$bin = as.numeric(as.character(metcov$bin))
metcov$start=c(1, head(metcov$bin+1,-1))
metcov = head(metcov,-1)

source('~/Tmp/AdditionalFile3_Code_DNase_seq.R')

gr.mii = GRanges(chr,IRanges(metcov$start, metcov$bin))
gr.mii$matrix = as.matrix(round(metcov[,grepl('MII', colnames(metcov)),with=FALSE],3))
AB.mii = createCorMatrix(gr.mii, method='spearman')
AB.mii = extractAB.dnase(AB.mii)

gr.mes = GRanges(chr,IRanges(metcov$start, metcov$bin))
gr.mes$matrix = as.matrix(round(metcov[,grepl('mES', colnames(metcov)),with=FALSE],3))
AB.mes = createCorMatrix(gr.mes, method='spearman')
AB.mes = extractAB.dnase(AB.mes)


library(AnnotationHub)
library(rtracklayer)
mii.mm10 = liftOver(AB.mii, chain = chain)
mii.mm10 = unlist(mii.mm10[ elementLengths(mii.mm10)==1])
mes.mm10 = liftOver(AB.mes, chain = chain)
mes.mm10 = unlist(mes.mm10[ elementLengths(mes.mm10)==1])



grmt = all_start
grmt$bin = as.numeric(as.character(grmt$bin))
grmt = grmt[order(grmt$bin),]

grmt = data.table(grmt)[,list(cnts = length(start)),by=list(bin, class)]
grmt = dcast(grmt, bin~class, value.var='cnts', fill=0)
grmt[,MT:=MT/sum(MT)]
grmt[,MT2:=MT2/sum(MT2)]

grmt$start=c(1, head(grmt$bin+1,-1))
grmt = na.omit(grmt)
ggrmt = GRanges(chr,IRanges(grmt$start, grmt$bin))
ggrmt$MT = grmt$MT
ggrmt$MT2 = grmt$MT2

fo = merge(data.table(as.matrix(findOverlaps(mii.mm10))),
           data.table(as.matrix(findOverlaps(mes.mm10))), by='queryHits')
fo$MT  = ggrmt$MT[fo$queryHits]
fo$MT2 = ggrmt$MT2[fo$queryHits]
fo$mii = mii.mm10$pc[fo$subjectHits.x]
fo$mes = mes.mm10$pc[fo$subjectHits.y]

fo = as.data.frame(fo[,-(1:3), with=F])
# fo[,1:2] = apply(fo[,1:2], 2, function(x)(x/(max(x))))
# fo[,3:4] = apply(fo[,3:4], 2, function(x)(x/(max(abs(x)))))
fo$rat = with(fo, log2((MT2+1)/(MT+1)))

foc = suppressWarnings(data.frame(id=1:nrow(fo),melt(fo)))

foc$type = foc$variable %in% c('mes','mii')
foc$type[foc$variable == 'rat'] = 'Ratio'

ggplot(data = foc, aes(x=id, y=value, color = variable)) + geom_density(stat='identity') + facet_wrap(~type, scales='free_y', ncol=1)

pdf('~/Tmp/plot8.pdf', width=42, height=12)
ggplot(data = foc, aes(x=id, y=value, color = variable)) + geom_density(stat='identity') + facet_wrap(~type, scales='free_y', ncol=1)
dev.off()