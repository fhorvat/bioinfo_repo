library(data.table)

load( "mirAnnot.dt.rda",verbose=T )

SHRINKLFC=TRUE
SHROF <- c(".ORIGINAL-LFC",".SHRUNKEN-LFC")[ as.integer(SHRINKLFC)+1 ]
if ( SHRINKLFC ) {
yLFCLIM=c(-7.5,+7.5)
} else           {
yLFCLIM=c(-9.5,+9.5)
}

#= EMB15.5
load(
 paste0("/storage/brno1-cerit/home/pepap/Taborska_smallRNAlibs_DicerXXE_191217/BAM/00.stats/01.miRNAcounts/01.final/EMB15.5-ES-smallRNAseq",SHROF,".de.dt.rda"),
 verbose=T
)

#= ESC
load(
 paste0("/storage/brno1-cerit/home/pepap/Eliska.Taborska/BAM/00.stats/01.miRNAcounts/01.final/ESC-ES-smallRNAseq",SHROF,".de.dt.rda"),
 verbose=T
)

#= TAR
load(
 paste0("/storage/brno1-cerit/home/pepap/Pullagura_2018_Genetics_PRJNA423238/BAM/READs/00.stats/01.miRNAcounts/01.final/EMB15.5-Pullagura-smallRNAseq.TAR",SHROF,".de.dt.rda"),
 verbose=T
)

#= GNT
load(
 paste0("/storage/brno1-cerit/home/pepap/Valeria.Buccheri/01.Dicer-helicase-domain-mutations/BAM.MERGED/00.stats/01.miRNAcounts/no.x19_KO/01.final/EMB15.5-DcrGNT-VB",SHROF,".de.dt.rda"),
 verbose=T
)

#>> EMB ~ ESC
plt.dt <-
 merge(
  x=get("EMB15.5-ES-smallRNAseq.de.dt")[,c("ID","log2FoldChange","padj","Name","miRNA.STR","mirtron","pep.type"),with=F],
  y=get("ESC-ES-smallRNAseq.de.dt")[    ,c("ID","log2FoldChange","padj")                                        ,with=F],
  by="ID",sort=F,all=T,suffixes=c(".emb",".esc")
 )
plt.dt                                                           <- merge( plt.dt,mirAnnot.pepType.dt[,c("ID","det.STR"),with=F],by="ID",sort=F,all=T )
plt.dt[["CEX"]]                                                  <- 1.0
plt.dt[ padj.emb<0.05 & padj.esc<0.05 ][["CEX"]]                 <- 1.5
plt.dt[["PCH"]]                                                  <- 16
plt.dt[ padj.emb<0.05 & padj.esc<0.05 & det.STR=="5p" ][["PCH"]] <- 25
plt.dt[ padj.emb<0.05 & padj.esc<0.05 & det.STR=="3p" ][["PCH"]] <- 24
plt.dt[["BG"]]                                                   <- "gray50"
plt.dt[ padj.emb<0.05 & padj.esc<0.05 & det.STR=="5p" ][["BG"]]  <- "blue"
plt.dt[ padj.emb<0.05 & padj.esc<0.05 & det.STR=="3p" ][["BG"]]  <- "red"
plt.dt[ padj.emb<0.05 & padj.esc<0.05 & mirtron==T    ][["BG"]]  <- "limegreen"
plt.dt[["COL"]]                                                  <- "gray50"
plt.dt[ padj.emb<0.05 & padj.esc<0.05 ][["COL"]]                 <- "transparent"
plt.dt                                                           <- plt.dt[ order(CEX) ]

xall <- cor.test( x=plt.dt[["log2FoldChange.emb"]],                                 y=plt.dt[["log2FoldChange.esc"]],                                 method="pearson" )
xsig <- cor.test( x=plt.dt[ padj.emb<0.05 & padj.esc<0.05 ][["log2FoldChange.emb"]],y=plt.dt[ padj.emb<0.05 & padj.esc<0.05 ][["log2FoldChange.esc"]],method="pearson" )

pdf( file=paste0("EMBvsESC--both-signif",SHROF,"-",pepDate(),".pdf"),width=25,height=25,pointsize=25 )
par( bg="white",pty="s" )
plot(
 plt.dt[,c("log2FoldChange.emb","log2FoldChange.esc"),with=F],
 col=plt.dt[["COL"]],cex=plt.dt[["CEX"]],pch=plt.dt[["PCH"]],bg=plt.dt[["BG"]],
 main="",panel.first=list(abline(v=0,h=0,col="gray50")),xlim=yLFCLIM,ylim=yLFCLIM,
 xlab="",
 ylab=""
)
mtext( text=expression( italic("Dicer"^{"X/X"})*"/wt E15.5" ),side=1,line=3,cex=1.5 )
mtext( text=expression( italic("Dicer"^{"X/X"})*"/wt ESC"   ),side=2,line=3,cex=1.5 )

text( x=yLFCLIM[1]*0.95,y=yLFCLIM[2]*0.95,labels=sprintf("R=%4.3f",xall$estimate),cex=1.1,col="black" )
text( x=yLFCLIM[1]*0.95,y=yLFCLIM[2]*0.90,labels=sprintf("R=%4.3f",xsig$estimate),cex=1.1,col="red"   )
dev.off()

plt.dt <-
 merge(
  x=get("EMB15.5-ES-smallRNAseq.de.dt")[,c("ID","log2FoldChange","padj","Name","miRNA.STR","mirtron","pep.type"),with=F],
  y=get("ESC-ES-smallRNAseq.de.dt")[    ,c("ID","log2FoldChange","padj")                                        ,with=F],
  by="ID",sort=F,all=T,suffixes=c(".emb",".esc")
 )
plt.dt                                                           <- merge( plt.dt,mirAnnot.pepType.dt[,c("ID","det.STR"),with=F],by="ID",sort=F,all=T )
plt.dt[["CEX"]]                                                  <- 1.0
plt.dt[ padj.emb<0.05 & padj.esc<1.05 ][["CEX"]]                 <- 1.5
plt.dt[["PCH"]]                                                  <- 16
plt.dt[ padj.emb<0.05 & padj.esc<1.05 & det.STR=="5p" ][["PCH"]] <- 25
plt.dt[ padj.emb<0.05 & padj.esc<1.05 & det.STR=="3p" ][["PCH"]] <- 24
plt.dt[["BG"]]                                                   <- "gray50"
plt.dt[ padj.emb<0.05 & padj.esc<1.05 & det.STR=="5p" ][["BG"]]  <- "blue"
plt.dt[ padj.emb<0.05 & padj.esc<1.05 & det.STR=="3p" ][["BG"]]  <- "red"
plt.dt[ padj.emb<0.05 & padj.esc<1.05 & mirtron==T    ][["BG"]]  <- "limegreen"
plt.dt[["COL"]]                                                  <- "gray50"
plt.dt[ padj.emb<0.05 & padj.esc<1.05 ][["COL"]]                 <- "transparent"
plt.dt                                                           <- plt.dt[ order(CEX) ]

xall <- cor.test( x=plt.dt[["log2FoldChange.emb"]],                                 y=plt.dt[["log2FoldChange.esc"]],                                 method="pearson" )
xsig <- cor.test( x=plt.dt[ padj.emb<0.05 & padj.esc<1.05 ][["log2FoldChange.emb"]],y=plt.dt[ padj.emb<0.05 & padj.esc<1.05 ][["log2FoldChange.esc"]],method="pearson" )

pdf( file=paste0("EMBvsESC--EMB-signif",SHROF,"-",pepDate(),".pdf"),width=25,height=25,pointsize=25 )
par( bg="white",pty="s" )
plot(
 plt.dt[,c("log2FoldChange.emb","log2FoldChange.esc"),with=F],
 col=plt.dt[["COL"]],cex=plt.dt[["CEX"]],pch=plt.dt[["PCH"]],bg=plt.dt[["BG"]],
 main="",panel.first=list(abline(v=0,h=0,col="gray50")),xlim=yLFCLIM,ylim=yLFCLIM,
 xlab="",
 ylab=""
)
mtext( text=expression( italic("Dicer"^{"X/X"})*"/wt E15.5" ),side=1,line=3,cex=1.5 )
mtext( text=expression( italic("Dicer"^{"X/X"})*"/wt ESC"   ),side=2,line=3,cex=1.5 )

text( x=yLFCLIM[1]*0.95,y=yLFCLIM[2]*0.95,labels=sprintf("R=%4.3f",xall$estimate),cex=1.1,col="black" )
text( x=yLFCLIM[1]*0.95,y=yLFCLIM[2]*0.90,labels=sprintf("R=%4.3f",xsig$estimate),cex=1.1,col="red"   )
dev.off()

#>> EMB ~ TAR
plt.dt <-
 merge(
  x=get("EMB15.5-ES-smallRNAseq.de.dt")[           ,c("ID","log2FoldChange","padj","Name","miRNA.STR","mirtron","pep.type"),with=F],
  y=get("EMB15.5-Pullagura-smallRNAseq.TAR.de.dt")[,c("ID","log2FoldChange","padj")                                        ,with=F],
  by="ID",sort=F,all=T,suffixes=c(".emb",".tar")
 )
plt.dt                                                           <- merge( plt.dt,mirAnnot.pepType.dt[,c("ID","det.STR"),with=F],by="ID",sort=F,all=T )
plt.dt[["CEX"]]                                                  <- 1.0
plt.dt[ padj.emb<1.05 & padj.tar<0.05 ][["CEX"]]                 <- 1.5
plt.dt[["PCH"]]                                                  <- 16
plt.dt[ padj.emb<1.05 & padj.tar<0.05 & det.STR=="5p" ][["PCH"]] <- 25
plt.dt[ padj.emb<1.05 & padj.tar<0.05 & det.STR=="3p" ][["PCH"]] <- 24
plt.dt[["BG"]]                                                   <- "gray50"
plt.dt[ padj.emb<1.05 & padj.tar<0.05 & det.STR=="5p" ][["BG"]]  <- "blue"
plt.dt[ padj.emb<1.05 & padj.tar<0.05 & det.STR=="3p" ][["BG"]]  <- "red"
plt.dt[ padj.emb<1.05 & padj.tar<0.05 & mirtron==T    ][["BG"]]  <- "limegreen"
plt.dt[["COL"]]                                                  <- "gray50"
plt.dt[ padj.emb<1.05 & padj.tar<0.05 ][["COL"]]                 <- "transparent"
plt.dt                                                           <- plt.dt[ order(CEX) ]

xall <- cor.test( x=plt.dt[["log2FoldChange.emb"]],                                 y=plt.dt[["log2FoldChange.tar"]],                                 method="pearson" )
xsig <- cor.test( x=plt.dt[ padj.emb<1.05 & padj.tar<0.05 ][["log2FoldChange.emb"]],y=plt.dt[ padj.emb<1.05 & padj.tar<0.05 ][["log2FoldChange.tar"]],method="pearson" )

pdf( file=paste0("EMBvsTAR--Tarbp2-signif",SHROF,"-",pepDate(),".pdf"),width=25,height=25,pointsize=25 )
par( bg="white",pty="s" )
plot(
 plt.dt[,c("log2FoldChange.emb","log2FoldChange.tar"),with=F],
 col=plt.dt[["COL"]],cex=plt.dt[["CEX"]],pch=plt.dt[["PCH"]],bg=plt.dt[["BG"]],
 main="",panel.first=list(abline(v=0,h=0,col="gray50")),xlim=yLFCLIM,ylim=yLFCLIM,
 xlab="",
 ylab=""
)
mtext( text=expression( italic("Dicer"^{"X/X"})*"/wt E15.5"  ),side=1,line=3,cex=1.5 )
mtext( text=expression( italic("Tarbp2"^{"-/-"})*"/wt E15.5" ),side=2,line=3,cex=1.5 )

text( x=yLFCLIM[1]*0.95,y=yLFCLIM[2]*0.95,labels=sprintf("R=%4.3f",xall$estimate),cex=1.1,col="black" )
text( x=yLFCLIM[1]*0.95,y=yLFCLIM[2]*0.90,labels=sprintf("R=%4.3f",xsig$estimate),cex=1.1,col="red"   )
dev.off()

#>> ESC ~ TAR
plt.dt <-
 merge(
  x=get("ESC-ES-smallRNAseq.de.dt")[               ,c("ID","log2FoldChange","padj","Name","miRNA.STR","mirtron","pep.type"),with=F],
  y=get("EMB15.5-Pullagura-smallRNAseq.TAR.de.dt")[,c("ID","log2FoldChange","padj")                                        ,with=F],
  by="ID",sort=F,all=T,suffixes=c(".esc",".tar")
 )
plt.dt                                                           <- merge( plt.dt,mirAnnot.pepType.dt[,c("ID","det.STR"),with=F],by="ID",sort=F,all=T )
plt.dt[["CEX"]]                                                  <- 1.0
plt.dt[ padj.esc<1.05 & padj.tar<0.05 ][["CEX"]]                 <- 1.5
plt.dt[["PCH"]]                                                  <- 16
plt.dt[ padj.esc<1.05 & padj.tar<0.05 & det.STR=="5p" ][["PCH"]] <- 25
plt.dt[ padj.esc<1.05 & padj.tar<0.05 & det.STR=="3p" ][["PCH"]] <- 24
plt.dt[["BG"]]                                                   <- "gray50"
plt.dt[ padj.esc<1.05 & padj.tar<0.05 & det.STR=="5p" ][["BG"]]  <- "blue"
plt.dt[ padj.esc<1.05 & padj.tar<0.05 & det.STR=="3p" ][["BG"]]  <- "red"
plt.dt[ padj.esc<1.05 & padj.tar<0.05 & mirtron==T    ][["BG"]]  <- "limegreen"
plt.dt[["COL"]]                                                  <- "gray50"
plt.dt[ padj.esc<1.05 & padj.tar<0.05 ][["COL"]]                 <- "transparent"
plt.dt                                                           <- plt.dt[ order(CEX) ]

xall <- cor.test( x=plt.dt[["log2FoldChange.esc"]],                                 y=plt.dt[["log2FoldChange.tar"]],                                 method="pearson" )
xsig <- cor.test( x=plt.dt[ padj.esc<1.05 & padj.tar<0.05 ][["log2FoldChange.esc"]],y=plt.dt[ padj.esc<1.05 & padj.tar<0.05 ][["log2FoldChange.tar"]],method="pearson" )

pdf( file=paste0("ESCvsTAR--Tarbp2-signif",SHROF,"-",pepDate(),".pdf"),width=25,height=25,pointsize=25 )
par( bg="white",pty="s" )
plot(
 plt.dt[,c("log2FoldChange.esc","log2FoldChange.tar"),with=F],
 col=plt.dt[["COL"]],cex=plt.dt[["CEX"]],pch=plt.dt[["PCH"]],bg=plt.dt[["BG"]],
 main="",panel.first=list(abline(v=0,h=0,col="gray50")),xlim=yLFCLIM,ylim=yLFCLIM,
 xlab="",
 ylab=""
)
mtext( text=expression( italic("Dicer"^{"X/X"})*"/wt ESC"    ),side=1,line=3,cex=1.5 )
mtext( text=expression( italic("Tarbp2"^{"-/-"})*"/wt E15.5" ),side=2,line=3,cex=1.5 )

text( x=yLFCLIM[1]*0.95,y=yLFCLIM[2]*0.95,labels=sprintf("R=%4.3f",xall$estimate),cex=1.1,col="black" )
text( x=yLFCLIM[1]*0.95,y=yLFCLIM[2]*0.90,labels=sprintf("R=%4.3f",xsig$estimate),cex=1.1,col="red"   )
dev.off()

#>> EMB ~ GNT
plt.dt <-
 merge(
  x=get("EMB15.5-ES-smallRNAseq.de.dt")[,c("ID","log2FoldChange","padj","Name","miRNA.STR","mirtron","pep.type"),with=F],
  y=get("EMB15.5-DcrGNT-VB.de.dt")[     ,c("ID","log2FoldChange","padj")                                        ,with=F],
  by="ID",sort=F,all=T,suffixes=c(".emb",".gnt")
 )
plt.dt                                                           <- merge( plt.dt,mirAnnot.pepType.dt[,c("ID","det.STR"),with=F],by="ID",sort=F,all=T )
plt.dt[["CEX"]]                                                  <- 1.0
plt.dt[ padj.emb<1.05 & padj.gnt<0.05 ][["CEX"]]                 <- 1.5
plt.dt[["PCH"]]                                                  <- 16
plt.dt[ padj.emb<1.05 & padj.gnt<0.05 & det.STR=="5p" ][["PCH"]] <- 25
plt.dt[ padj.emb<1.05 & padj.gnt<0.05 & det.STR=="3p" ][["PCH"]] <- 24
plt.dt[["BG"]]                                                   <- "gray50"
plt.dt[ padj.emb<1.05 & padj.gnt<0.05 & det.STR=="5p" ][["BG"]]  <- "blue"
plt.dt[ padj.emb<1.05 & padj.gnt<0.05 & det.STR=="3p" ][["BG"]]  <- "red"
plt.dt[ padj.emb<1.05 & padj.gnt<0.05 & mirtron==T    ][["BG"]]  <- "limegreen"
plt.dt[["COL"]]                                                  <- "gray50"
plt.dt[ padj.emb<1.05 & padj.gnt<0.05 ][["COL"]]                 <- "transparent"
plt.dt                                                           <- plt.dt[ order(CEX) ]

xall <- cor.test( x=plt.dt[["log2FoldChange.emb"]],                                 y=plt.dt[["log2FoldChange.gnt"]],                                 method="pearson" )
xsig <- cor.test( x=plt.dt[ padj.emb<1.05 & padj.gnt<0.05 ][["log2FoldChange.emb"]],y=plt.dt[ padj.emb<1.05 & padj.gnt<0.05 ][["log2FoldChange.gnt"]],method="pearson" )

pdf( file=paste0("EMBvsGNT--only-GNT-signif",SHROF,"-",pepDate(),".pdf"),width=25,height=25,pointsize=25 )
par( bg="white",pty="s" )
plot(
 plt.dt[,c("log2FoldChange.emb","log2FoldChange.gnt"),with=F],
 col=plt.dt[["COL"]],cex=plt.dt[["CEX"]],pch=plt.dt[["PCH"]],bg=plt.dt[["BG"]],
 main="",panel.first=list(abline(v=0,h=0,col="gray50")),xlim=yLFCLIM,ylim=yLFCLIM,
 xlab="",
 ylab=""
)
mtext( text=expression( italic("Dicer"^{"X/X"})*"/wt E15.5"     ),side=1,line=3,cex=1.5 )
mtext( text=expression( italic("Dicer"^{"GNT/GNT"})*"/wt E15.5" ),side=2,line=3,cex=1.5 )

text( x=yLFCLIM[1]*0.95,y=yLFCLIM[2]*0.95,labels=sprintf("R=%4.3f",xall$estimate),cex=1.1,col="black" )
text( x=yLFCLIM[1]*0.95,y=yLFCLIM[2]*0.90,labels=sprintf("R=%4.3f",xsig$estimate),cex=1.1,col="red"   )
dev.off()

