library(GenomicAlignments)
library(rtracklayer)
library(data.table)

#>>> miR-3102; miR-15a; miR-361; miR-145a; miR-323; miR-291b
preLoc <-
 GRanges(
  "chr7:100882307-100882409","chr14:61632038-61632096","chrX:113074827-113074888","chr18:61647827-61647888","chr12:109712523-109712578","chr7:3219494-3219550"
  preID = c( "MI0014099","MI0000564","MI0000761","MI0000169","MI0000592","MI0003539" )
 )
PREIDS <- mcols(preLoc)[["preID"]]
mirLoc <- 
 GRanges(
  c(
   "chr7:100882388-100882409","chr7:100882307-100882329",
   "chr14:61632075-61632096","chr14:61632038-61632059","chrX:113074867-113074888","chrX:113074827-113074851","chr18:61647866-61647888","chr18:61647827-61647848",
   "chr12:109712523-109712544","chr12:109712558-109712578","chr7:3219494-3219515","chr7:3219529-3219550"
  ),
  Name  = c(
   "mmu-miR-3102-5p","mmu-miR-3102-3p",
   "mmu-miR-15a-5p","mmu-miR-15a-3p","mmu-miR-361-5p","mmu-miR-361-3p","mmu-miR-145a-5p","mmu-miR-145a-3p",
   "mmu-miR-323-5p","mmu-miR-323-3p","mmu-miR-291b-5p","mmu-miR-291b-3p"
   ),
  ID    = c(
   "MIMAT0014933","MIMAT0014936",
   "MIMAT0000526","MIMAT0004624","MIMAT0000704","MIMAT0017075","MIMAT0000157","MIMAT0004534",
   "MIMAT0004638","MIMAT0000551","MIMAT0003189","MIMAT0003190"
   ),
  preID = c(
   "MI0014099","MI0014099",
   "MI0000564","MI0000564","MI0000761","MI0000761","MI0000169","MI0000169",
   "MI0000592","MI0000592","MI0003539","MI0003539"
   )
 )

rngRL<-c(21,23)
rngNH<-c(1,100)
rngHI<-c(1,100)
rngNM<-c(0,0)
minOVR=15
minRN=10

BAMS <-
 c(
  "RS10_Mos_r1.bam",
  "RS10_Mos_r2.bam",
  "RS10_Mos_r3.bam",
  "RS7_Mos_r1.bam",
  "RS7_Mos_r2.bam",
  "RS7_Mos_r3.bam"
 )
NAME <- gsub("[.]bam","",BAMS)

#>>> FUNCTIONS
plot.one <-
 function(
  inCV,add.annot=F,yMAX=NULL,xLAB=T,yLAB=T,xTIT=NULL,xLWD=10,preCOL="black",mirCOL="gray50",xNFC=1,
  preGR=preLoc,mirGR=mirLoc,CONST=1e06
 ) {

 preCV <- inCV[ preGR ][[1]]
 mirGR <- mirGR[ order(start(mirGR)) ]

 if (xLAB) {xLAB=as.character(preGR)} else {xLAB=""}
 if (yLAB) {yLAB="Expression (RPM)"}     else {yLAB=""}

 plot.dt <-
  data.table(
   pos=seq(start(preGR),end(preGR),1),
   as.data.table( preCV ),
   col=preCOL
  )
 plot.dt[ pos %in% unlist(sapply(X=seq(length(mirGR)),FUN=function(x){ return( seq(start(mirGR[x]),end(mirGR[x]),1) ) })) ][["col"]] <- mirCOL
 plot.dt[["value"]] <- plot.dt[["value"]]*xNFC
 if (is.null(yMAX)) { yMAX=max(plot.dt[,value])*1.0 }
 if (add.annot) { yMIN=(-3)*(yMAX/50) } else { yMIN=0 }
 plot(
  plot.dt[,c("pos","value"),with=F], type="h", xlab=xLAB, ylab=yLAB, main=xTIT, font.lab=2, xlim=range(plot.dt[,pos])+c(-10,10), ylim=c(yMIN,yMAX), bty="n", lend="square",lwd=xLWD, col=plot.dt[,col]
 )
 if (add.annot) {
  arrows( x0=start(preGR), x1=end(preGR), y0=yMIN*0.5, y1=yMIN*0.5, lwd=05, col=preCOL, code=0 )
  rect( xleft=start(mirGR)-0.5, xright=end(mirGR)+0.5, ybottom=yMIN*0.7, ytop=yMIN*0.3, col=mirCOL )
  text( x=abs(end(mirGR)+start(mirGR))/2, y=yMIN*1.1,labels=mcols(mirGR)[,"Name"] )
 }

 return(plot.dt)

}
#<<<

i               <- 1
total.mapped.dt <- data.table()
for ( iB in BAMS ) {

 cat("\n ++ ",iB,"\n",sep="")
 tmp.ga                                      <- readGAlignments(file=iB,param=ScanBamParam(what="qname",tag=c("NH","HI","nM")))
 tmp.gr                                      <- as(tmp.ga,"GRanges")
 tmp.gr                                      <- tmp.gr[
  ( width(tmp.gr) >= rngRL[1] )  & ( width(tmp.gr) <= rngRL[2] )  &
  mcols(tmp.gr)[,"NH"]>=rngNH[1] & mcols(tmp.gr)[,"NH"]<=rngNH[2] &
  mcols(tmp.gr)[,"HI"]>=rngHI[1] & mcols(tmp.gr)[,"HI"]<=rngHI[2] &
  mcols(tmp.gr)[,"nM"]>=rngNM[1] & mcols(tmp.gr)[,"nM"]<=rngNM[2]
 ]
 total.mapped.dt                             <- rbind( total.mapped.dt,data.table(NAME=NAME[i],total=sum( !duplicated(mcols(tmp.gr)[,"qname"]) )) )
 tmp.ol                                      <- findOverlaps( query=tmp.gr,subject=mature_miRNAs_mmu.gr,minoverlap=minOVR,ignore.strand=F)
 tmp.dt                                      <-
  data.table(
   ID=mcols(mature_miRNAs_mmu.gr)[subjectHits(tmp.ol),"ID"],
   NH=mcols(              tmp.gr)[  queryHits(tmp.ol),"NH"]
  )
 tmp.dt                                      <- tmp.dt[,{ list(cnt=length(NH),frc=sum(1/NH)) },by="ID"]
 colnames(tmp.dt)[ colnames(tmp.dt)=="cnt" ] <- paste(NAME[i],".cnt",sep="")
 colnames(tmp.dt)[ colnames(tmp.dt)=="frc" ] <- paste(NAME[i],".frc",sep="")
 tmp.gr                                      <- tmp.gr[ unique(queryHits(tmp.ol)) ]

 c <- 1
 for ( inh in sort(unique(mcols(tmp.gr)[,"NH"])) ) {
  if ( c==1 ) {
   tmp.cv <-   coverage( tmp.gr[ mcols(tmp.gr)[,"NH"]==inh ] )*(1/inh)
  } else                   {
   tmp.cv <- ( coverage( tmp.gr[ mcols(tmp.gr)[,"NH"]==inh ] )*(1/inh) ) + tmp.cv
  }
  c <- (c+1)
 }

 assign(x=paste(NAME[i],".gr",sep=""),value=tmp.gr)
 assign(x=paste(NAME[i],".cv",sep=""),value=tmp.cv)

 rm( list="tmp.cv" )

 i <- (i+1)

}

ODIR1="FIGs_mean/"
ODIR2="FIGs_mean_noScale/"
ODIR3="FIGs_mean_combined/"
dir.create(ODIR1)
dir.create(ODIR2)
dir.create(ODIR3)
CONST=1e06
rs07cv01 <- (get(paste(NAME[05],".cv",sep=""))/total.mapped.dt[ NAME==NAME[05],total ])*CONST
rs07cv02 <- (get(paste(NAME[06],".cv",sep=""))/total.mapped.dt[ NAME==NAME[06],total ])*CONST
rs07cv03 <- (get(paste(NAME[07],".cv",sep=""))/total.mapped.dt[ NAME==NAME[07],total ])*CONST
rs10cv01 <- (get(paste(NAME[01],".cv",sep=""))/total.mapped.dt[ NAME==NAME[01],total ])*CONST
rs10cv02 <- (get(paste(NAME[02],".cv",sep=""))/total.mapped.dt[ NAME==NAME[02],total ])*CONST
rs10cv03 <- (get(paste(NAME[03],".cv",sep=""))/total.mapped.dt[ NAME==NAME[03],total ])*CONST

for ( preid in PREIDS ) {

curr.ymax <- max(c( (( rs07cv01+rs07cv02+rs07cv03 )/3)[ preLoc[ mcols(preLoc)[["preID"]]==preid ] ][[1]],(( rs10cv01+rs10cv02+rs10cv03 )/3)[ preLoc[ mcols(preLoc)[["preID"]]==preid ] ][[1]] ))

pdf(file=paste(ODIR1,preid,".pdf",sep=""),width=40,height=60,pointsize=50)
 par(bg="white",mfrow=c(2,1))
 plot.one(
  inCV = ( rs07cv01+rs07cv02+rs07cv03 )/3,add.annot=T,yMAX=curr.ymax,xLAB=T,xTIT="RS7 : mean coverage", xNFC=1,xLWD=25,
  preGR=preLoc[ mcols(preLoc)[["preID"]]==preid ],mirGR=mirLoc[ mcols(mirLoc)[["preID"]]==preid ]
 )
 plot.one(
  inCV = ( rs10cv01+rs10cv02+rs10cv03 )/3,add.annot=T,yMAX=curr.ymax,xLAB=T,xTIT="RS10 : mean coverage",xNFC=1,xLWD=25,
  preGR=preLoc[ mcols(preLoc)[["preID"]]==preid ],mirGR=mirLoc[ mcols(mirLoc)[["preID"]]==preid ]
 )
dev.off()

pdf(file=paste(ODIR2,preid,".pdf",sep=""),width=40,height=60,pointsize=50)
 par(bg="white",mfrow=c(2,1))
 plot.one(
  inCV = ( rs07cv01+rs07cv02+rs07cv03 )/3,add.annot=T,yMAX=NULL,     xLAB=T,xTIT="RS7 : mean coverage", xNFC=1,xLWD=25,
  preGR=preLoc[ mcols(preLoc)[["preID"]]==preid ],mirGR=mirLoc[ mcols(mirLoc)[["preID"]]==preid ]
 )
 plot.one(
  inCV = ( rs10cv01+rs10cv02+rs10cv03 )/3,add.annot=T,yMAX=NULL,     xLAB=T,xTIT="RS10 : mean coverage",xNFC=1,xLWD=25,
  preGR=preLoc[ mcols(preLoc)[["preID"]]==preid ],mirGR=mirLoc[ mcols(mirLoc)[["preID"]]==preid ]
 )
dev.off()

pdf(file=paste(ODIR3,preid,".pdf",sep=""),width=80,height=60,pointsize=50)
 par(bg="white",mfrow=c(2,2))
 plot.one(
  inCV = ( rs07cv01+rs07cv02+rs07cv03 )/3,add.annot=T,yMAX=curr.ymax,xLAB=T,xTIT="RS7 : mean coverage", xNFC=1,xLWD=25,
  preGR=preLoc[ mcols(preLoc)[["preID"]]==preid ],mirGR=mirLoc[ mcols(mirLoc)[["preID"]]==preid ]
 )
 plot.one(
  inCV = ( rs07cv01+rs07cv02+rs07cv03 )/3,add.annot=T,yMAX=NULL,     xLAB=T,xTIT="RS7 : mean coverage", xNFC=1,xLWD=25,
  preGR=preLoc[ mcols(preLoc)[["preID"]]==preid ],mirGR=mirLoc[ mcols(mirLoc)[["preID"]]==preid ]
 )
 plot.one(
  inCV = ( rs10cv01+rs10cv02+rs10cv03 )/3,add.annot=T,yMAX=curr.ymax,xLAB=T,xTIT="RS10 : mean coverage",xNFC=1,xLWD=25,
  preGR=preLoc[ mcols(preLoc)[["preID"]]==preid ],mirGR=mirLoc[ mcols(mirLoc)[["preID"]]==preid ]
 )
 plot.one(
  inCV = ( rs10cv01+rs10cv02+rs10cv03 )/3,add.annot=T,yMAX=NULL,     xLAB=T,xTIT="RS10 : mean coverage",xNFC=1,xLWD=25,
  preGR=preLoc[ mcols(preLoc)[["preID"]]==preid ],mirGR=mirLoc[ mcols(mirLoc)[["preID"]]==preid ]
 )
dev.off()

}

