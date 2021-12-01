# change if necessary:
nthreads=10

library(data.table)
setDTthreads(nthreads)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(stringdist)
library(stringr)

#Testing testing:
#Note that for now this works on 

getsequences <- function(gr, windowsize=50, LEN=150, SOMEVALUE=55){
  seqstart <- start(gr)-5000
  chrs <- as.character(seqnames(gr))
  
  length <- end(gr)+5000-seqstart
  starts <- seq(seqstart, seqstart+length,windowsize)
  ends <- starts+LEN
  
  x <- GRanges(chrs, IRanges(starts, ends))
  somesequence <- getSeq(BSgenome.Mmusculus.UCSC.mm10,x)
  alldists<-as.data.table(stringdistmatrix(somesequence, somesequence))
  alldists[,row:=starts]
  alldists <- melt(alldists, id.vars=c("row"))
  alldists[,':='(col=ends[as.numeric(str_extract(variable,"\\d+"))], chr=unique(chrs))]
  alldists[value==0, value:=1000]
  alldists <- alldists[row<col,.(chr,row,col, value)][order(row)]
  
  x <- alldists[value<SOMEVALUE]
  x[,follows:={k=(row==dplyr::lag(row)+windowsize);
  l=(col==dplyr::lag(col)+windowsize);
  p=ifelse(k==l & l==1,1,0)
  p = rleid(p)
  ifelse(p%%2==0,p,(p+1))%/%2
  }]
  x[,.(chr=unique(chr),rowstart=min(row),rowend=max(row),colstart=min(col),colend=max(col),value=mean(value)),follows]
}


lll <- fread("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/Documentation/L1s_nested_ours_20180516.ZJM.tidy.csv")
gr <-  GRanges(lll)[1]
sek <- getsequences(gr)
grseq <- GRanges(sek$chr,IRanges(sek$rowstart, sek$colend))
countOverlaps(GRanges(sek$chr,IRanges(sek$rowstart, sek$colend)), gr, ignore.strand=T)


