library(data.table)
library(GenomicRanges)
library(DESeq2)
library(pcaExplorer)
library(matrixStats)
library(PTXQC)

load( "premirAnnot.dt.rda", verbose=T )

OFNAME="EMB15.5-ES-smallRNAseq"
eMIN=100
SHRINKLFC=TRUE
SHROF <- c(".ORIGINAL-LFC",".SHRUNKEN-LFC")[ as.integer(SHRINKLFC)+1 ]
if ( SHRINKLFC ) {
yLFCLIM=c(-7.5,+7.5)
} else           {
yLFCLIM=c(-9.5,+9.5)
}
xEXPLIM=c(0.0,6.0)

RM_COLS=c("ID","Chr","Start","End","Strand","Length")

load( "mirAnnot.dt.rda",verbose=T )
mirAnnot.dt[["type"]] <- mirAnnot.dt[["pep.type"]]

#>>> FUNCTIONS
featureCountTab2dt <- function(fcTableFile,libname.prefix="Rep_") {

 cat(" ++ Loading counting table from featureCount\n",sep="")
 if ( substr( x = fcTableFile, start = nchar(fcTableFile)-1, stop = nchar(fcTableFile) ) == "gz" ) {
  out.dt <- fread( cmd = paste("zcat ",fcTableFile," | grep -v \"^#\"",sep=""), header = T, stringsAsFactors = F )
#  out.dt <- fread( input = paste("zcat ",fcTableFile," | grep -v \"^#\"",sep=""), header = T, stringsAsFactors = F )
 } else                                                                                            {
  out.dt <- fread( cmd = paste("cat ", fcTableFile," | grep -v \"^#\"",sep=""), header = T, stringsAsFactors = F )
#  out.dt <- fread( input = paste("cat ", fcTableFile," | grep -v \"^#\"",sep=""), header = T, stringsAsFactors = F )
 }

 cat(" ++ Checking expected format of the table\n",sep="")
 cat("     Expected column names : \"Geneid\" \"Chr\" \"Start\" \"End\" \"Strand\" \"Length\"\n",sep="")
 cat("     YourFile column names : \"",paste(colnames(out.dt)[1:6],collapse="\" \""),"\"\n",sep="")
 if ( any( colnames(out.dt)[1:6] != c("Geneid","Chr","Start","End","Strand","Length") ) ) {
  warning(" Changes in the featureCount output format ??\n Need to be checked manually !!")
 } else                                                      {
  colnames(out.dt)[1] <- "GeneID"
 }

 cat(" ++ Basic formating of the column names\n",sep="")
 colnames(out.dt)[ seq(7,length(colnames(out.dt))) ] <-
  gsub(
   pattern     = "^.*/",
   replacement = "",
   x           = colnames(out.dt)[ seq(7,length(colnames(out.dt))) ]
  )

 #= Remove all (at least 4-character) strings in BAM file names
 while( nchar(LCSn( strings = colnames(out.dt)[ seq(7,length(colnames(out.dt))) ], min_LCS_length = 4 )) >= 4 ) {
  colnames(out.dt)[ seq(7,length(colnames(out.dt))) ] <-
   gsub(
    pattern     = LCSn( strings = colnames(out.dt)[ seq(7,length(colnames(out.dt))) ], min_LCS_length = 4 ),
    replacement = "",
    x           = colnames(out.dt)[ seq(7,length(colnames(out.dt))) ]
   )
 }

 colnames(out.dt)[ seq(7,length(colnames(out.dt))) ] <-
  gsub( pattern = "[._]$", replacement = "", x = colnames(out.dt)[ seq(7,length(colnames(out.dt))) ] )

 colnames(out.dt)[ seq(7,length(colnames(out.dt))) ] <-
  paste(libname.prefix,colnames(out.dt)[ seq(7,length(colnames(out.dt))) ],sep="")

 cat( " ++ Finished\n",sep="")
 return(out.dt)

}
mrgByPreMiRNA <- function( idt,LIMS=c(-6.5,+6.5),SELE=sele.miRNAs ) {

 odt <- merge.data.table( premirAnnot.dt,idt[,c("ID","log2FoldChange","padj"),with=F],by.x="ID.5p",by.y="ID",all.x=T,sort=F                         )
 odt <- merge.data.table(            odt,idt[,c("ID","log2FoldChange","padj"),with=F],by.x="ID.3p",by.y="ID",all.x=T,sort=T,suffixes=c(".5p",".3p") )

 odt <- odt[ ( ( type.5p=="miRNA" ) & ( type.3p=="miRNA*" ) ) | ( ( type.5p=="miRNA*" ) & ( type.3p=="miRNA" ) ) ]

 odt <- odt[ ( ID.5p %in% SELE ) | ( ID.3p %in% SELE ) ]

 odt[["PCH"]] <- 16
 odt[ type.5p=="miRNA" ][["PCH"]] <- 25
 odt[ type.3p=="miRNA" ][["PCH"]] <- 24
 odt[["COL"]] <- "gray50"
 odt[ type.5p=="miRNA" ][["COL"]] <- "blue"
 odt[ type.3p=="miRNA" ][["COL"]] <- "red"
 odt[ type.5p=="miRNA" & ( padj.5p>=0.05 & padj.3p>=0.05 ) ][["COL"]] <- rgb(0,0,100,alpha=15,maxColorValue=100)
 odt[ type.3p=="miRNA" & ( padj.5p>=0.05 & padj.3p>=0.05 ) ][["COL"]] <- rgb(100,0,0,alpha=15,maxColorValue=100)
 odt <- odt[ order(COL,decreasing=F) ]
 odt[["CEX"]] <- 1.75

 odt[["X"]] <- as.numeric("1000")
 odt[ type.5p=="miRNA"  ][["X"]] <- odt[ type.5p=="miRNA"  ][["log2FoldChange.5p"]]
 odt[ type.3p=="miRNA"  ][["X"]] <- odt[ type.3p=="miRNA"  ][["log2FoldChange.3p"]]
 odt[["Y"]] <- as.numeric("1000")
 odt[ type.5p=="miRNA*" ][["Y"]] <- odt[ type.5p=="miRNA*" ][["log2FoldChange.5p"]]
 odt[ type.3p=="miRNA*" ][["Y"]] <- odt[ type.3p=="miRNA*" ][["log2FoldChange.3p"]]

 par( bg="white",pty="s" )
 plot(
  odt[,c("X","Y"),with=F],
  panel.first=abline(h=0,v=0,col="gray75",lty=1,lwd=1),font.lab=2,
  xlab="miRNA",ylab="miRNA*",
  xlim=LIMS,ylim=LIMS,
  col=rgb(0,0,0,alpha=25,maxColorValue=100),
  pch=odt[["PCH"]],bg =odt[["COL"]],cex=odt[["CEX"]]
 )

 return( odt )

}
pep.mirMAplot <- function( MIRDT,xTYPE="miRNA",TIT=NULL,lfcLIM=NULL,xLIM=NULL,mirtron=F,COLS=c("blue","red","gray50","limegreen") ) {

 cat( "\n miRNA type : \"",xTYPE,"\"\n",sep="")
 cat( " Acceptable inputs : \"miRNA\", \"miRNA*\", \"all\"\n",sep="" )

 if ( is.null(TIT) ) {
  TIT <- sub( "[*]","-passenger",x=TIT )
 }

if ( xTYPE=="all" ) {
 MIRDT          <- MIRDT[ type %in% unique(type) ]
} else {
 MIRDT          <- MIRDT[ type==xTYPE ]
}
 MIRDT[["LOG"]] <- log10(MIRDT[["baseMean"]]+1)
 print( range(MIRDT[["LOG"]]) )
 print( range(MIRDT[["log2FoldChange"]]) )
 MIRDT[["BG"]]                    <- COLS[3]
 MIRDT[ miRNA.STR=="5p" ][["BG"]] <- COLS[1]
 MIRDT[ miRNA.STR=="3p" ][["BG"]] <- COLS[2]
 if ( mirtron ) {
 MIRDT[ mirtron==T ][["BG"]] <- COLS[4]
 }

 if ( is.null(lfcLIM) ) {
 YMAX           <- max(abs(na.omit(MIRDT[["log2FoldChange"]])))
 YLIM           <- c(YMAX*(-1.01),YMAX*(+1.01))
 } else                 {
 YLIM           <- lfcLIM
 }

 if ( is.null(xLIM) ) {
 XLIM           <- range(MIRDT[["LOG"]])*c(0.99,1.01)
 } else               {
 XLIM           <- xLIM
 }
 nonsig.PCH=16
 nonsig.CEX=1.50
 nonsig.COL=COLS[3]
 nonsig.BG ="transparent"
 par( bg="white" )
 plot(
  MIRDT[,c("LOG","log2FoldChange"),with=F],
  col=nonsig.COL,bg=nonsig.BG,pch=nonsig.PCH,cex=nonsig.CEX,
  xlab="mean of normalized counts (log10)",ylab="log2FoldChange",main=TIT,xlim=XLIM,ylim=YLIM,
  xaxt="n",
  font.lab=2
 )
 axis( side=1,at=seq(XLIM[1],XLIM[2],1),labels=c( expression(" 1"),expression("10"),expression("10"^"2"),expression("10"^"3"),expression("10"^"4"),expression("10"^"5"),expression("10"^"6") ) )
 sig5.PCH=25
 sig5.CEX=2.25
 sig5.COL="transparent"
 sig5.BG =COLS[1]
 points(
  MIRDT[ miRNA.STR=="5p" & padj<=0.05 ,c("LOG","log2FoldChange"),with=F],
  pch=sig5.PCH,cex=sig5.CEX,col=sig5.COL,bg=MIRDT[ miRNA.STR=="5p" & padj<=0.05 , BG ]
 )
 sig3.PCH=24
 sig3.CEX=2.25
 sig3.COL="transparent"
 sig3.BG =COLS[2]
 points(
  MIRDT[ miRNA.STR=="3p" & padj<=0.05 ,c("LOG","log2FoldChange"),with=F],
  pch=sig3.PCH,cex=sig3.CEX,col=sig3.COL,bg=MIRDT[ miRNA.STR=="3p" & padj<=0.05 , BG ]
 )
 legend( "topleft",pch=c(sig5.PCH,sig3.PCH),legend=c("5'","3'"),bty="n",col=c(sig5.COL,sig3.COL),pt.bg=c(sig5.BG,sig3.BG),pt.cex=c(sig5.CEX,sig3.CEX),y.intersp=2 )

 return( MIRDT )

}
#<<<

fc_all.dt <- featureCountTab2dt( fcTableFile=FCTABLEFILE,libname.prefix="" )

xREF="WT"

exp.mat <- as.matrix( x=fc_all.dt[ , c("ID",colnames(fc_all.dt)[ !( colnames(fc_all.dt) %in% RM_COLS ) ]) , with=F ],rownames="ID" )
exp.df  <- data.frame( CON=as.factor( sub("_.*$","",sub("^DicerX_","",colnames(exp.mat))) ),WHO=as.factor("ES"),row.names=colnames(exp.mat) )

FC.dds      <-
 DESeqDataSetFromMatrix(
  countData  = round( exp.mat ),
  colData    = exp.df,
  design     = ~ CON
 )
FC.dds$CON  <- relevel( x=FC.dds$CON,ref=xREF )
FC.dds      <- DESeq( FC.dds )

KO2WT.dt <- dds2res( DDS=FC.dds,CON="CON",A="KO",B=xREF,as.DT=T,GIDCOL="ID",shrink=SHRINKLFC,shrink.type="normal" )
KO2WT.dt <- merge( KO2WT.dt,mirAnnot.dt,by="ID",all.x=T,sort=T )

EMB.exp.dt <- merge( mirAnnot.pepType.dt[,c("ID","miRNA.STR","det.STR","mirtron","Name","pep.type"),with=F],as.data.table(counts( FC.dds,normalized=T ),keep.rownames="ID"),by="ID",sort=F,all=T )
EMB.exp.dt <- merge( EMB.exp.dt,KO2WT.dt[,c("ID","log2FoldChange","padj"),with=F],by="ID",sort=F,all=T )
colnames(EMB.exp.dt)[ colnames(EMB.exp.dt)=="log2FoldChange" ] <- "log2FoldChange.EMB"
colnames(EMB.exp.dt)[ colnames(EMB.exp.dt)=="padj"           ] <- "padj.EMB"
save( EMB.exp.dt,file="EMB.exp.dt.rda" )

print( yLFCLIM )
print( range(KO2WT.dt[["log2FoldChange"]]) )
if (T) {

if (T) {
pdf( file=paste(OFNAME,SHROF,".DESeq2.MAplot-",pepDate(),".pdf",sep=""),         width=50,height=20,pointsize=25 )
par( bg="white",mfrow=c(1,2) )
pep.mirMAplot( MIRDT=KO2WT.dt,mirtron=F,xTYPE="miRNA", lfcLIM=yLFCLIM,xLIM=xEXPLIM,TIT=paste0(OFNAME,"-main")      )
pep.mirMAplot( MIRDT=KO2WT.dt,mirtron=F,xTYPE="miRNA*",lfcLIM=yLFCLIM,xLIM=xEXPLIM,TIT=paste0(OFNAME,"-passenger") )
dev.off()

pdf( file=paste(OFNAME,SHROF,".DESeq2.MAplot-mirtron-",pepDate(),".pdf",sep=""), width=50,height=20,pointsize=25 )
par( bg="white",mfrow=c(1,2) )
pep.mirMAplot( MIRDT=KO2WT.dt,mirtron=T,xTYPE="miRNA", lfcLIM=yLFCLIM,xLIM=xEXPLIM,TIT=paste0(OFNAME,"-main")      )
pep.mirMAplot( MIRDT=KO2WT.dt,mirtron=T,xTYPE="miRNA*",lfcLIM=yLFCLIM,xLIM=xEXPLIM,TIT=paste0(OFNAME,"-passenger") )
dev.off()
}

pdf( file=paste(OFNAME,SHROF,".DESeq2.MAplot-all-",pepDate(),".pdf",sep=""),     width=25,height=20,pointsize=25 )
par( bg="white",mfrow=c(1,1) )
pep.mirMAplot( MIRDT=KO2WT.dt,mirtron=T,xTYPE="all",   lfcLIM=yLFCLIM,xLIM=xEXPLIM,TIT=paste0(OFNAME,"-all")      )
dev.off()

ddd <- KO2WT.dt[ ,c("ID","Name","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","pep.type","mirtron"),with=F ]
ddd <- ddd[ order(padj) ]
rrr <- merge( as.data.table(fpm( object=FC.dds,robust=F ),keep.rownames="ID"),KO2WT.dt[,c("ID","padj","log2FoldChange","Name","pep.type","mirtron"),with=F],by="ID",sort=F,all=T )
rrr <- rrr[ order(padj) ]
save( list=c("ddd","rrr"),file=paste0("CSVs",SHROF,".rda") )

assign( x=paste0(OFNAME,".de.dt"),value=KO2WT.dt )
save( list=paste0(OFNAME,".de.dt"),file=paste0(OFNAME,SHROF,".de.dt.rda") )

}

if (T) {

ALL.dds      <-
 DESeqDataSetFromMatrix(
  countData  = round( exp.mat ),
  colData    = exp.df,
  design     = ~ 1
 )

xPCA.df <- pcaExplorer::pcaplot( varianceStabilizingTransformation(ALL.dds),        intgroup=c("CON"),     text_labels=T, title=OFNAME, point_size=10.0, returnData=T )
save( xPCA.df,file="xPCA.df.rda" )

xCEX=3.5
pdf( file=paste(OFNAME,".DESeq2.PCA-",pepDate(),".pdf",sep=""),                  width=25,height=20,pointsize=25 )
par( bg="white",mar=c(5,5,5,5),xpd=T )
plot(
 xPCA.df[,c("PC1","PC2")],pch=16,col=c("gray50","red")[ as.integer(xPCA.df$CON=="KO")+1 ],cex=xCEX,
 xlab="PC1 (67.0% variance)",ylab="PC2 (13.6% variance)",main=OFNAME
)
legend( x=max(xPCA.df$PC1)*1.1,y=max(xPCA.df$PC2)*0.9,legend=c("WT\n","KO"),bty="n",pch=16,pt.cex=xCEX,col=c("gray50","red") )
dev.off()

}

load( paste0("all.DE",SHROF,".dt.rda"),verbose=T )
sele.miRNAs <- unique(all.DE.dt[ ( baseMean.EMB>=eMIN ) | ( baseMean.ESC>=eMIN ) , ID ])

pdf( file=paste0(OFNAME,SHROF,"-crosshair-",pepDate(),".pdf"),                   width=25,height=25,pointsize=25 )
qqq <- mrgByPreMiRNA( idt=KO2WT.dt,LIMS=yLFCLIM,SELE=sele.miRNAs )
dev.off()

