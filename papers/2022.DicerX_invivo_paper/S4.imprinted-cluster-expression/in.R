library(data.table)
library(GenomicRanges)

load( "mirAnnot.dt.rda",verbose=T )

#= EMB15.5-ES
emb15.5.DE.dt <- load( "/storage/brno1-cerit/home/pepap/Taborska_smallRNAlibs_DicerXXE_191217/BAM/00.stats/01.miRNAcounts/EMB15.5-ES-smallRNAseq.de.dt.rda",verbose=T )
emb15.5.DE.dt <- get(emb15.5.DE.dt)

#= ESC-ES
esc.DE.dt <- load( "/storage/brno1-cerit/home/pepap/Eliska.Taborska/BAM/00.stats/01.miRNAcounts/ESC-ES-smallRNAseq.de.dt.rda",verbose=T )
esc.DE.dt <- get(esc.DE.dt)

#= EMB15.5-Pu
pkr.DE.dt <- load( "/storage/brno1-cerit/home/pepap/Pullagura_2018_Genetics_PRJNA423238/BAM/READs/00.stats/01.miRNAcounts/EMB15.5-Pullagura-smallRNAseq.PKR.de.dt.rda",verbose=T )
pkr.DE.dt <- get(pkr.DE.dt)
tar.DE.dt <- load( "/storage/brno1-cerit/home/pepap/Pullagura_2018_Genetics_PRJNA423238/BAM/READs/00.stats/01.miRNAcounts/EMB15.5-Pullagura-smallRNAseq.TAR.de.dt.rda",verbose=T )
tar.DE.dt <- get(tar.DE.dt)

#= EMB15.5-DcrGNT-VB
gnt.DE.dt <- load( "/storage/brno1-cerit/home/pepap/Valeria.Buccheri/01.Dicer-helicase-domain-mutations/BAM.MERGED/00.stats/01.miRNAcounts/no.x19_KO/EMB15.5-DcrGNT-VB.de.dt.rda",verbose=T )
gnt.DE.dt <- get(gnt.DE.dt)

#= FUN
getChrReg <-
 function(
  chReg=c("chr12:109704935-109749758"),
  deDT=c("emb15.5.DE.dt","esc.DE.dt","pkr.DE.dt",  "tar.DE.dt",  "gnt.DE.dt"),
  deNM=c("EMB15.5-ES",   "ESC-ES",   "EMB15.5-PKR","EMB15.5-TAR","DcrGNT-VB"),
  oPCH=c(24,25,21),oSTR=c("3p","5p","unknown"),
  oCEX=c(0.50,1.00,1.50,2.00),
  oxBY=10000
 ) {

 xCHR <- as.character( seqnames(GRanges(chReg)) )
 xBEG <-                  start(GRanges(chReg))
 xEND <-                    end(GRanges(chReg))

 xMIR <- mirAnnot.dt[ ( chr==xCHR ) & ( start >= xBEG ) & ( end <= xEND ) ][ order(start) ]

 par( bg="white",mar=c(3,10,3,1) )
 plot(
  x=-100,y=-100,xlim=c(0,nrow(xMIR)+1),ylim=c(0,length(deDT)+1),
  xaxt="n",yaxt="n",bty="n",
  xlab="",ylab="",main=chReg
 )
 axis( side=2,at=seq(length(deDT)),labels=rev(deNM),las=2,font=2 )
 xTIC <- seq( from=xBEG,to=xEND,by=oxBY )
 axis( side=1,at=seq( from=1,to=nrow(xMIR),length.out=length(xTIC) ),labels=round(xTIC),las=1,font=1 )

 j=0
 for ( idt in deDT ) {

  idt  <- merge( xMIR[,c("ID"),with=F],get(idt),by="ID",sort=F,all.x=T )

  xPOS <- seq_along(idt[["start"]])

  yPOS <- rep.int( (length(deDT)-j),nrow(idt) )

  xPCH                                <- rep.int( oPCH[3],nrow(idt) )
  xPCH[ idt[["miRNA.STR"]]==oSTR[1] ] <- oPCH[1]
  xPCH[ idt[["miRNA.STR"]]==oSTR[2] ] <- oPCH[2]

  xCEX <- c(NA,oCEX[1])[ as.integer(   idt[["baseMean"]] <= 10                                     )+1 ]
  xCEX[                              ( idt[["baseMean"]] >  10   ) & ( idt[["baseMean"]] <= 100  )     ] <- oCEX[2]
  xCEX[                              ( idt[["baseMean"]] >  100  ) & ( idt[["baseMean"]] <= 1000 )     ] <- oCEX[3]
  xCEX[                              ( idt[["baseMean"]] >  1000 )                                     ] <- oCEX[4]

  xBG  <- c("blue","red")[ as.integer( idt[["log2FoldChange"]] >= 0 )+1 ]
  xBG[ idt[["padj"]] >  0.05 ] <- "gray50"

  xPDJ <- c("","*")[ as.integer( idt[["padj"]] <= 0.05 )+1 ]
   
  points( x=xPOS,y=yPOS,    pch=xPCH,cex=xCEX,bg=xBG )
  points( x=xPOS,y=yPOS+0.2,pch=xPDJ,cex=1.00        )

  j <- j+1

 }

 return( xMIR )

}
#=

pdf( file=paste0("chr12.miRNAs-",pepDate(),".pdf"),width=40,height=15,pointsize=20 )
 getChrReg( chReg=c("chr12:109704935-109749758") )
dev.off()

