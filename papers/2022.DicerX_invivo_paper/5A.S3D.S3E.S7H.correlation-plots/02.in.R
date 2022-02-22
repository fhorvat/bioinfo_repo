library(data.table)

load( "mirAnnot.dt.rda",verbose=T )

SHRINKLFC=TRUE
SHROF <- c(".ORIGINAL-LFC",".SHRUNKEN-LFC")[ as.integer(SHRINKLFC)+1 ]
if ( SHRINKLFC ) {
yLFCLIM=c(-7.5,+7.5)
} else           {
yLFCLIM=c(-9.5,+9.5)
}
yLFCLIM=c(0,6)

#= EMB15.5
load(
 paste0("/storage/brno1-cerit/home/pepap/Taborska_smallRNAlibs_DicerXXE_191217/BAM/00.stats/01.miRNAcounts/01.final/CSVs",SHROF,".rda"),
 verbose=T
)
plt.dt             <- rrr[,c("ID","DicerX_WT_1_r3.frac","DicerX_WT_6_r1.frac","DicerX_WT_8B_r2.frac","padj","Name","pep.type","mirtron"),with=F]
colnames(plt.dt)[ colnames(plt.dt)=="padj" ] <- "padj.EMB"
plt.dt[["EMB.wt"]] <- rowMeans( as.matrix(plt.dt[,c("DicerX_WT_1_r3.frac","DicerX_WT_6_r1.frac","DicerX_WT_8B_r2.frac"),with=F]) )

#= ESC
load(
 paste0("/storage/brno1-cerit/home/pepap/Eliska.Taborska/BAM/00.stats/01.miRNAcounts/01.final/CSVs",SHROF,".rda"),
 verbose=T
)
plt.dt             <- merge( plt.dt,rrr[,c("ID","s_RS7_Mos_r1.frac","s_RS7_Mos_r2.frac","s_RS7_Mos_r3.frac","padj"),with=F] )
colnames(plt.dt)[ colnames(plt.dt)=="padj" ] <- "padj.ESC"
plt.dt[["ESC.wt"]] <- rowMeans( as.matrix(plt.dt[,c("s_RS7_Mos_r1.frac","s_RS7_Mos_r2.frac","s_RS7_Mos_r3.frac"),with=F]) )

#= Pullagura
load(
 paste0("/storage/brno1-cerit/home/pepap/Pullagura_2018_Genetics_PRJNA423238/BAM/READs/00.stats/01.miRNAcounts/01.final/CSVs",SHROF,".rda"),
 verbose=T
)
plt.dt             <- merge( plt.dt,rrr1[,c("ID","embryo_Tarbp2_WT_r1.frac","embryo_Tarbp2_WT_r2.frac","embryo_Tarbp2_WT_r3.frac","padj"),with=F] )
colnames(plt.dt)[ colnames(plt.dt)=="padj" ] <- "padj.TAR"
plt.dt[["TAR.wt"]] <- rowMeans( as.matrix(plt.dt[,c("embryo_Tarbp2_WT_r1.frac","embryo_Tarbp2_WT_r2.frac","embryo_Tarbp2_WT_r3.frac"),with=F]) )
plt.dt             <- merge( plt.dt,rrr2[,c("ID","embryo_Prkra_WT_r1.frac","embryo_Prkra_WT_r2.frac","padj"),with=F] )
colnames(plt.dt)[ colnames(plt.dt)=="padj" ] <- "padj.PKR"
plt.dt[["PKR.wt"]] <- rowMeans( as.matrix(plt.dt[,c("embryo_Prkra_WT_r1.frac","embryo_Prkra_WT_r2.frac"),with=F]) )

#= GNT
load(
 paste0("/storage/brno1-cerit/home/pepap/Valeria.Buccheri/01.Dicer-helicase-domain-mutations/BAM.MERGED/00.stats/01.miRNAcounts/no.x19_KO/01.final/CSVs",SHROF,".rda"),
 verbose=T
)
plt.dt             <- merge( plt.dt,rrr[,c("ID","x11_WT.frac","x14_WT.frac","x16_WT.frac","padj"),with=F] )
colnames(plt.dt)[ colnames(plt.dt)=="padj" ] <- "padj.GNT"
plt.dt[["GNT.wt"]] <- rowMeans( as.matrix(plt.dt[,c("x11_WT.frac","x14_WT.frac","x16_WT.frac"),with=F]) )

plt.dt             <- merge( plt.dt,mirAnnot.pepType.dt[,c("ID","det.STR"),with=F],by="ID",sort=F,all=T )
plt.RPM.dt         <- plt.dt
save( plt.RPM.dt,file="plt.RPM.dt.rda" )

#= FUN

getPTS <- function( SETS=c("EMB","ESC"),PLIM=c(0.05,0.05),LOG10=T,INPDT=plt.dt ) {

 IDS          <-        INPDT[ ( get(paste0("padj.",SETS[1]))<PLIM[1] ) & ( get(paste0("padj.",SETS[2]))<PLIM[2] ) ][["ID"]]
 STR          <-        INPDT[ ( get(paste0("padj.",SETS[1]))<PLIM[1] ) & ( get(paste0("padj.",SETS[2]))<PLIM[2] ) ][["det.STR"]]
 MRT          <-        INPDT[ ( get(paste0("padj.",SETS[1]))<PLIM[1] ) & ( get(paste0("padj.",SETS[2]))<PLIM[2] ) ][["mirtron"]]
 XXX          <- log10( INPDT[ ( get(paste0("padj.",SETS[1]))<PLIM[1] ) & ( get(paste0("padj.",SETS[2]))<PLIM[2] ) ][[paste0(SETS[1],".wt")]]+1 )
 YYY          <- log10( INPDT[ ( get(paste0("padj.",SETS[1]))<PLIM[1] ) & ( get(paste0("padj.",SETS[2]))<PLIM[2] ) ][[paste0(SETS[2],".wt")]]+1 )
 COR          <- cor.test( x=XXX,y=YYY,method="pearson" )

 PCH          <- c(   24,    25)[ as.integer(STR=="5p")+1 ]
 BG           <- c("red","blue")[ as.integer(STR=="5p")+1 ]
 BG[ MRT==T ] <- "limegreen"

 OUT          <- list( A=XXX,B=YYY,C=COR$estimate,D=COR$p.value,E=PCH,FF=BG )
 names(OUT)   <- c( SETS,"pearson","p.value","PCH","BG" )

 return( OUT )

}

#=

#>> EMB ~ ESC

CONS <- setNames( object=c("E15.5","ESC","Tarbp2","Prkra","GNT"),nm=c("EMB","ESC","TAR","PKR","GNT") )
ILST <- list( c("EMB","ESC"),c("EMB","ESC") )
IPVL <- list( c( 0.05,0.05 ),c( 0.05, 1.05) )

for ( i in seq_along(ILST) ) {

 print( ILST[[i]] )

PLOTDT <- getPTS( SETS=ILST[[i]],PLIM=IPVL[[i]],LOG10=T,INPDT=plt.dt )
SIGNIF <- c( c("-nonsig","-signif")[ as.integer(IPVL[[i]][1]==0.05)+1 ],c("-nonsig","-signif")[ as.integer(IPVL[[i]][2]==0.05)+1 ] )
pdf( file=paste0("wtExp--RPM--",ILST[[i]][1],SIGNIF[1],"-vs-",ILST[[i]][2],SIGNIF[2],SHROF,"-",pepDate(),".pdf"),width=25,height=25,pointsize=25 )
par( bg="white",pty="s" )
plot(
 x=PLOTDT[[ ILST[[i]][1] ]],y=PLOTDT[[ ILST[[i]][2] ]],
 col="transparent",bg=PLOTDT[["BG"]],cex=1.5,pch=PLOTDT[["PCH"]],
 main=paste0(ILST[[i]][1],SIGNIF[1],"-vs-",ILST[[i]][2],SIGNIF[2]),
 panel.first=list(abline(b=1,a=0,col="black",lty=2)),xlim=yLFCLIM,ylim=yLFCLIM,
 xlab="",
 ylab=""
)
XL <- CONS[ ILST[[i]][1] ]
YL <- CONS[ ILST[[i]][2] ]
mtext( text=bquote( "log"[10]*"(wt"~.(XL)*")" ),side=1,line=3,cex=1.5 )
mtext( text=bquote( "log"[10]*"(wt"~.(YL)*")" ),side=2,line=3,cex=1.5 )

text( x=yLFCLIM[1]+0.2,y=yLFCLIM[2]*0.95,labels=sprintf("R=%4.3f",PLOTDT[["pearson"]]),cex=1.1,col="black" )
dev.off()

}

