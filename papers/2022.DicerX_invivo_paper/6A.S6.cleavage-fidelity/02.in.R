library(data.table)
library(matrixStats)
library(pheatmap)
library(openxlsx)

load( "DicerX-CPs.19to25nt-20220109.rda", verbose=T )
load( "mirAnnot-CPs.dt.rda",              verbose=T )
load( "all.DE.dt.rda",                    verbose=T )

eMIN=100
ONAME="19to25nt"
TOPN=50

breakList <- seq(-1,+1,0.01)

EMB.mirnas <- all.DE.dt[ baseMean.EMB > eMIN , ID ]
ESC.mirnas <- all.DE.dt[ baseMean.ESC > eMIN , ID ]
TAR.mirnas <- all.DE.dt[ baseMean.TAR > eMIN , ID ]
PKR.mirnas <- all.DE.dt[ baseMean.PKR > eMIN , ID ]
GNT.mirnas <- all.DE.dt[ baseMean.GNT > eMIN , ID ]

EMB.ES.CON <- c(
  "DicerX_KO_10B_r2","DicerX_HET_12B_r2","DicerX_HET_13B_r3","DicerX_WT_1_r3","DicerX_HET_2B_r4","DicerX_KO_2_r5","DicerX_KO_3B_r3",
  "DicerX_KO_4_r1","DicerX_WT_6_r1","DicerX_KO_7B_r4","DicerX_HET_7_r1","DicerX_WT_8B_r2"
 )
ESC.ES.CON <- c(
 "RS10_Mos_r1","RS10_Mos_r2","RS10_Mos_r3","RS10_NC_r1",
 "RS7_Mos_r1", "RS7_Mos_r2", "RS7_Mos_r3", "RS7_NC_r1",
 "RSP_Mos_r1", "RSP_Mos_r2", "RSP_Mos_r3", "RSP_NC_r1"
 )
EMB.PU.CON <- c(
 "emb_Prkra_Mut_r1", "emb_Prkra_Mut_r2", "emb_Prkra_Mut_r3",
 "emb_Prkra_WT_r1",  "emb_Prkra_WT_r2",
 "emb_Tarbp2_Mut_r1","emb_Tarbp2_Mut_r2","emb_Tarbp2_Mut_r3",
 "emb_Tarbp2_WT_r1", "emb_Tarbp2_WT_r2", "emb_Tarbp2_WT_r3"
 )
GNT.VB.CON <- c(
 "x11_WT","x14_WT","x16_WT",
 "x18_GNT_heterozygot","x19_GNT_homozygot",
 "x21","x22","x23","x24",
 "x3_GNT_homozygot","x4_GNT_homozygot","x9_GNT_homozygot"
 )

MATNAMES <- c( EMB.ES.CON,ESC.ES.CON,EMB.PU.CON,GNT.VB.CON )

#= FCE
# (1) matrix lines => density
normMat <- function( iMAT ) {

 out.mat <- apply(
  X      = iMAT,
  MARGIN = 1,
  FUN    = function(L){
   if ( sum(L)==0 ) {
    return(L)
   } else           {
    return(L/sum(L))
   }
  }
  )

 return( t(out.mat) )

}

# (2) define new CP based on the expression; report multiple CPs & OUTlier CPs
deffCP <- function( iMAT,SEL=as.character(seq(-10,+10)),OUT=5,ANNOT=mirAnnot.dt ) {

 all.rwn <- rownames( iMAT )
 ADD     <- as.integer(SEL[1]) - as.integer(colnames(iMAT)[1])
 CANON   <- seq_along(SEL)[ SEL=="0" ]

 j=1
 out.mat <-
  apply(
   X      = iMAT,
   MARGIN = 1,
   FUN    = function(L){
    LSEL <- L[ SEL ]
    tmp.rnw <- rownames(iMAT)[j]
    j <<- ( j + 1 )
    if ( sum(LSEL)!=0 ) {
     ZERO <- max( LSEL )
     ZERO <- seq_along(LSEL)[ LSEL==ZERO ]
     if ( length(ZERO)>1 ) {
      print( c(ZERO,tmp.rnw) )
      return( CANON+ADD )
     }
     if ( abs(CANON-ZERO)>OUT ) {
      print( c(ZERO,tmp.rnw) )
      return( CANON+ADD )
     } else               {
      return( ZERO +ADD )
     }
    } else              {
     return(  CANON+ADD )
    }
   }
  )

 return( out.mat[all.rwn] )

}

# (3) rename colnames in matrix, center "0" to new CP
deffCPcenter <- function( iMAT,CP.DT,EXPAND=5 ) {

 out.mat <-
  sapply(
   X   = seq(nrow(CP.DT)),
   FUN = function(i){
    cLINE        <- iMAT[ CP.DT[ i , ID ] , seq( from=CP.DT[ i , CP-EXPAND ],to=CP.DT[ i , CP+EXPAND ] ) ]
    names(cLINE) <- as.character(seq((-1)*EXPAND,(+1)*EXPAND))
    return( cLINE )
   },
   simplify=T,
   USE.NAMES=T
  )

 out.mat           <- t(out.mat)
 rownames(out.mat) <- CP.DT[ , ID ]

 return( out.mat )

}

# (4) for each genomic position return "FCE"-value of ( As - Bs )
diffMat <- function( As,Bs,FCE="median" ) {

 i=1
 tmp.ncol <- c()
 for ( a in As ) {
 for ( b in Bs ) {
  if ( i == 1 )  {
   tmp.mat <- ( get(a) - get(b) )
  } else         {
   tmp.mat <- cbind( tmp.mat,( get(a) - get(b) ) )
  }
  tmp.ncol <- append( tmp.ncol,ncol(get(a)) )
  tmp.ncol <- append( tmp.ncol,ncol(get(b)) )
  i <- ( i + 1 )
 }
 }

 tmp.ncol <- unique( tmp.ncol )
# print( tmp.ncol )
 if ( length(tmp.ncol)!=1 ) { stop("\n! Input matrices differ in number of columns !\n\n") }

 rNAME    <- rownames( tmp.mat )
 cNAME    <- c()
 out.vals <- c()
 for ( i in c( seq(1,(tmp.ncol-1)),0 ) ) {
  cNAME    <- append( cNAME,unique(colnames(tmp.mat)[ ( seq_along(tmp.mat[1,]) %% tmp.ncol ) == i ]) )
  out.vals <- append( out.vals,apply( X=tmp.mat,MARGIN=1,FUN=function(L){ return(get(FCE)( L[ ( seq_along(L) %% tmp.ncol ) == i ] )) } ) )
 }

 out.mat <- matrix( data=out.vals,byrow=F,ncol=tmp.ncol,dimnames=list(rNAME,cNAME) )

 return( out.mat )

}

#=

if ( T ) {

#= I. define miRNA cleavage-points (CPs) from WT reads
#= (1) EMB-ES
 EMB.CPs.dt <- data.table(
  ID=EMB.mirnas
 )
 EMB.CPs.dt <- merge( EMB.CPs.dt,mirAnnot.dt[,c("ID","gLoc","preID.gLoc","Name","det.STR","type"),with=F],by="ID",all.x=T,sort=F )
 EMB.CPs.dt <- merge( EMB.CPs.dt,  all.DE.dt[,c("ID","baseMean.EMB"),                             with=F],by="ID",all.x=T,sort=F )

 tmp.dt <- deffCP( iMAT=DicerX_WT_1_r3.CPs.mat[EMB.mirnas,], SEL=as.character(seq(-10,+10)),OUT=7,ANNOT=mirAnnot.dt )
 tmp.dt <- data.table( ID=names(tmp.dt),DicerX_WT_1_r3=tmp.dt )
 EMB.CPs.dt <- merge( EMB.CPs.dt,tmp.dt,by="ID",sort=F,all=T )

 tmp.dt <- deffCP( iMAT=DicerX_WT_6_r1.CPs.mat[EMB.mirnas,], SEL=as.character(seq(-10,+10)),OUT=7,ANNOT=mirAnnot.dt )
 tmp.dt <- data.table( ID=names(tmp.dt),DicerX_WT_6_r1=tmp.dt )
 EMB.CPs.dt <- merge( EMB.CPs.dt,tmp.dt,by="ID",sort=F,all=T )
#= ?? manually checked : "MIMAT0001342"

 tmp.dt <- deffCP( iMAT=DicerX_WT_8B_r2.CPs.mat[EMB.mirnas,],SEL=as.character(seq(-10,+10)),OUT=7,ANNOT=mirAnnot.dt )
 tmp.dt <- data.table( ID=names(tmp.dt),DicerX_WT_8B_r2=tmp.dt )
 EMB.CPs.dt <- merge( EMB.CPs.dt,tmp.dt,by="ID",sort=F,all=T )
#= ?? manually checked : "MIMAT0003508"

 EMB.CPs.dt[["CP"]] <- rowMedians(as.matrix(EMB.CPs.dt[,c("DicerX_WT_1_r3","DicerX_WT_6_r1","DicerX_WT_8B_r2"),with=F]))

#= (2) ESC-ES
 ESC.CPs.dt <- data.table(
  ID=ESC.mirnas
 )
 ESC.CPs.dt <- merge( ESC.CPs.dt,mirAnnot.dt[,c("ID","gLoc","preID.gLoc","Name","det.STR","type"),with=F],by="ID",all.x=T,sort=F )
 ESC.CPs.dt <- merge( ESC.CPs.dt,  all.DE.dt[,c("ID","baseMean.ESC"),                             with=F],by="ID",all.x=T,sort=F )

 tmp.dt <- deffCP( iMAT=RS7_Mos_r1.CPs.mat[ESC.mirnas,],SEL=as.character(seq(-10,+10)),OUT=7,ANNOT=mirAnnot.dt )
 tmp.dt <- data.table( ID=names(tmp.dt),RS7_Mos_r1=tmp.dt )
 ESC.CPs.dt <- merge( ESC.CPs.dt,tmp.dt,by="ID",sort=F,all=T )

 tmp.dt <- deffCP( iMAT=RS7_Mos_r2.CPs.mat[ESC.mirnas,],SEL=as.character(seq(-10,+10)),OUT=7,ANNOT=mirAnnot.dt )
 tmp.dt <- data.table( ID=names(tmp.dt),RS7_Mos_r2=tmp.dt )
 ESC.CPs.dt <- merge( ESC.CPs.dt,tmp.dt,by="ID",sort=F,all=T )

 tmp.dt <- deffCP( iMAT=RS7_Mos_r3.CPs.mat[ESC.mirnas,],SEL=as.character(seq(-10,+10)),OUT=7,ANNOT=mirAnnot.dt )
 tmp.dt <- data.table( ID=names(tmp.dt),RS7_Mos_r3=tmp.dt )
 ESC.CPs.dt <- merge( ESC.CPs.dt,tmp.dt,by="ID",sort=F,all=T )
#= ?? manually checked : "MIMAT0001081", "MIMAT0031410"

 ESC.CPs.dt[["CP"]] <- rowMedians(as.matrix(ESC.CPs.dt[,c("RS7_Mos_r1","RS7_Mos_r2","RS7_Mos_r3"),with=F]))

#= (3) PKR-Pu
 PKR.CPs.dt <- data.table(
  ID=PKR.mirnas
 )
 PKR.CPs.dt <- merge( PKR.CPs.dt,mirAnnot.dt[,c("ID","gLoc","preID.gLoc","Name","det.STR","type"),with=F],by="ID",all.x=T,sort=F )
 PKR.CPs.dt <- merge( PKR.CPs.dt,  all.DE.dt[,c("ID","baseMean.PKR"),                             with=F],by="ID",all.x=T,sort=F )

 tmp.dt <- deffCP( iMAT=emb_Prkra_WT_r1.CPs.mat[PKR.mirnas,],SEL=as.character(seq(-10,+10)),OUT=7,ANNOT=mirAnnot.dt )
 tmp.dt <- data.table( ID=names(tmp.dt),emb_Prkra_WT_r1=tmp.dt )
 PKR.CPs.dt <- merge( PKR.CPs.dt,tmp.dt,by="ID",sort=F,all=T )

 tmp.dt <- deffCP( iMAT=emb_Prkra_WT_r2.CPs.mat[PKR.mirnas,],SEL=as.character(seq(-10,+10)),OUT=7,ANNOT=mirAnnot.dt )
 tmp.dt <- data.table( ID=names(tmp.dt),emb_Prkra_WT_r2=tmp.dt )
 PKR.CPs.dt <- merge( PKR.CPs.dt,tmp.dt,by="ID",sort=F,all=T )
#= ?? manually checked : "MIMAT0017005", "MIMAT0003473", "MIMAT0003473_1"

 PKR.CPs.dt[["CP"]] <- rowMedians(as.matrix(PKR.CPs.dt[,c("emb_Prkra_WT_r1","emb_Prkra_WT_r2"),with=F]))

#= (4) TAR-Pu
 TAR.CPs.dt <- data.table(
  ID=TAR.mirnas
 )
 TAR.CPs.dt <- merge( TAR.CPs.dt,mirAnnot.dt[,c("ID","gLoc","preID.gLoc","Name","det.STR","type"),with=F],by="ID",all.x=T,sort=F )
 TAR.CPs.dt <- merge( TAR.CPs.dt,  all.DE.dt[,c("ID","baseMean.TAR"),                             with=F],by="ID",all.x=T,sort=F )

 tmp.dt <- deffCP( iMAT=emb_Tarbp2_WT_r1.CPs.mat[TAR.mirnas,],SEL=as.character(seq(-10,+10)),OUT=7,ANNOT=mirAnnot.dt )
 tmp.dt <- data.table( ID=names(tmp.dt),emb_Tarbp2_WT_r1=tmp.dt )
 TAR.CPs.dt <- merge( TAR.CPs.dt,tmp.dt,by="ID",sort=F,all=T )

 tmp.dt <- deffCP( iMAT=emb_Tarbp2_WT_r2.CPs.mat[TAR.mirnas,],SEL=as.character(seq(-10,+10)),OUT=7,ANNOT=mirAnnot.dt )
 tmp.dt <- data.table( ID=names(tmp.dt),emb_Tarbp2_WT_r2=tmp.dt )
 TAR.CPs.dt <- merge( TAR.CPs.dt,tmp.dt,by="ID",sort=F,all=T )

 tmp.dt <- deffCP( iMAT=emb_Tarbp2_WT_r3.CPs.mat[TAR.mirnas,],SEL=as.character(seq(-10,+10)),OUT=7,ANNOT=mirAnnot.dt )
 tmp.dt <- data.table( ID=names(tmp.dt),emb_Tarbp2_WT_r3=tmp.dt )
 TAR.CPs.dt <- merge( TAR.CPs.dt,tmp.dt,by="ID",sort=F,all=T )
#= ?? manually checked : "MIMAT0016982"

 TAR.CPs.dt[["CP"]] <- rowMedians(as.matrix(TAR.CPs.dt[,c("emb_Tarbp2_WT_r1","emb_Tarbp2_WT_r2","emb_Tarbp2_WT_r3"),with=F]))

#= (5) GNT-VB
 GNT.CPs.dt <- data.table(
  ID=GNT.mirnas
 )
 GNT.CPs.dt <- merge( GNT.CPs.dt,mirAnnot.dt[,c("ID","gLoc","preID.gLoc","Name","det.STR","type"),with=F],by="ID",all.x=T,sort=F )
 GNT.CPs.dt <- merge( GNT.CPs.dt,  all.DE.dt[,c("ID","baseMean.GNT"),                             with=F],by="ID",all.x=T,sort=F )

 tmp.dt <- deffCP( iMAT=x11_WT.CPs.mat[GNT.mirnas,],SEL=as.character(seq(-10,+10)),OUT=7,ANNOT=mirAnnot.dt )
 tmp.dt <- data.table( ID=names(tmp.dt),x11_WT=tmp.dt )
 GNT.CPs.dt <- merge( GNT.CPs.dt,tmp.dt,by="ID",sort=F,all=T )

 tmp.dt <- deffCP( iMAT=x14_WT.CPs.mat[GNT.mirnas,],SEL=as.character(seq(-10,+10)),OUT=7,ANNOT=mirAnnot.dt )
 tmp.dt <- data.table( ID=names(tmp.dt),x14_WT=tmp.dt )
 GNT.CPs.dt <- merge( GNT.CPs.dt,tmp.dt,by="ID",sort=F,all=T )

 tmp.dt <- deffCP( iMAT=x16_WT.CPs.mat[GNT.mirnas,],SEL=as.character(seq(-10,+10)),OUT=7,ANNOT=mirAnnot.dt )
 tmp.dt <- data.table( ID=names(tmp.dt),x16_WT=tmp.dt )
 GNT.CPs.dt <- merge( GNT.CPs.dt,tmp.dt,by="ID",sort=F,all=T )
#= ?? manually checked : "MIMAT0014805"

 GNT.CPs.dt[["CP"]] <- rowMedians(as.matrix(GNT.CPs.dt[,c("x11_WT","x14_WT","x16_WT"),with=F]))

#= II. center CP-matrices to new CPs & normalize
 out.objs    <- c()
 for ( iM in paste(EMB.ES.CON,".CPs.mat",sep="") ) {

  cat( " => ",iM,"\n",sep="" )
  tmp.mat <- normMat( iMAT=deffCPcenter( iMAT=get(iM),CP.DT=EMB.CPs.dt,EXPAND=5 ) )
  print(           dim(tmp.mat)  )
  print( table(rowSums(tmp.mat)) )
  assign( x=sub( "[.]CPs[.]mat$",".CPs.normCent.mat",iM ),value=tmp.mat )
  out.objs <- append( out.objs,sub( "[.]CPs[.]mat$",".CPs.normCent.mat",iM ) )
 
 }
 for ( iM in paste(ESC.ES.CON,".CPs.mat",sep="") ) {

  cat( " => ",iM,"\n",sep="" )
  tmp.mat <- normMat( iMAT=deffCPcenter( iMAT=get(iM),CP.DT=ESC.CPs.dt,EXPAND=5 ) )
  print(           dim(tmp.mat)  )
  print( table(rowSums(tmp.mat)) )
  assign( x=sub( "[.]CPs[.]mat$",".CPs.normCent.mat",iM ),value=tmp.mat )
  out.objs <- append( out.objs,sub( "[.]CPs[.]mat$",".CPs.normCent.mat",iM ) )
 
 }
 for ( iM in paste(EMB.PU.CON[1:5],".CPs.mat",sep="") ) {

  cat( " => ",iM,"\n",sep="" )
  tmp.mat <- normMat( iMAT=deffCPcenter( iMAT=get(iM),CP.DT=PKR.CPs.dt,EXPAND=5 ) )
  print(           dim(tmp.mat)  )
  print( table(rowSums(tmp.mat)) )
  assign( x=sub( "[.]CPs[.]mat$",".CPs.normCent.mat",iM ),value=tmp.mat )
  out.objs <- append( out.objs,sub( "[.]CPs[.]mat$",".CPs.normCent.mat",iM ) )
 
 }
 for ( iM in paste(EMB.PU.CON[6:11],".CPs.mat",sep="") ) {

  cat( " => ",iM,"\n",sep="" )
  tmp.mat <- normMat( iMAT=deffCPcenter( iMAT=get(iM),CP.DT=TAR.CPs.dt,EXPAND=5 ) )
  print(           dim(tmp.mat)  )
  print( table(rowSums(tmp.mat)) )
  assign( x=sub( "[.]CPs[.]mat$",".CPs.normCent.mat",iM ),value=tmp.mat )
  out.objs <- append( out.objs,sub( "[.]CPs[.]mat$",".CPs.normCent.mat",iM ) )
 
 }
 for ( iM in paste(GNT.VB.CON,".CPs.mat",sep="") ) {

  cat( " => ",iM,"\n",sep="" )
  tmp.mat <- normMat( iMAT=deffCPcenter( iMAT=get(iM),CP.DT=GNT.CPs.dt,EXPAND=5 ) )
  print(           dim(tmp.mat)  )
  print( table(rowSums(tmp.mat)) )
  assign( x=sub( "[.]CPs[.]mat$",".CPs.normCent.mat",iM ),value=tmp.mat )
  out.objs <- append( out.objs,sub( "[.]CPs[.]mat$",".CPs.normCent.mat",iM ) )
 
 }
 save( list=out.objs,file=paste("DicerX-CPs.normCent.",ONAME,"-",pepDate(),".rda",sep="") )
# load( file="DicerX-CPs.normCent.19to25nt-20220110.rda",verbose=T )

#= III. unify KO/WT differences as a median-value for each genomic position
 EMB.diff <-
  diffMat(
   As=paste(c("DicerX_KO_10B_r2","DicerX_KO_2_r5","DicerX_KO_3B_r3","DicerX_KO_4_r1","DicerX_KO_7B_r4"),".CPs.normCent.mat",sep=""),
   Bs=paste(c("DicerX_WT_1_r3",  "DicerX_WT_6_r1","DicerX_WT_8B_r2"),                                   ".CPs.normCent.mat",sep=""),
   FCE="median"
  )
# >> order from -1 => +1 at "0" position ( == increasing )
 EMB.diff   <- EMB.diff[   order(   (EMB.diff[,"0"]),   decreasing=F) , ]
 tmp.dt     <- data.table( EMB.diff,keep.rownames="ID" )
 EMB.CPs.dt <- merge( EMB.CPs.dt,tmp.dt,by="ID",all=T )
# >> order from biggest changes (abs.values) to smallest at "0" position ( == decreasing )
 EMB.CPs.dt <- EMB.CPs.dt[ order(abs(EMB.CPs.dt[["0"]]),decreasing=T) ]

 ESC.diff <-
  diffMat(
   As=paste(c("RS10_Mos_r1","RS10_Mos_r2","RS10_Mos_r3"),".CPs.normCent.mat",sep=""),
   Bs=paste(c("RS7_Mos_r1", "RS7_Mos_r2", "RS7_Mos_r3"), ".CPs.normCent.mat",sep=""),
   FCE="median"
  )
# >> order from -1 => +1 at "0" position ( == increasing )
 ESC.diff   <- ESC.diff[   order(   (ESC.diff[,"0"]),   decreasing=F) , ]
 tmp.dt     <- data.table( ESC.diff,keep.rownames="ID" )
 ESC.CPs.dt <- merge( ESC.CPs.dt,tmp.dt,by="ID",all=T )
# >> order from biggest changes (abs.values) to smallest at "0" position ( == decreasing )
 ESC.CPs.dt <- ESC.CPs.dt[ order(abs(ESC.CPs.dt[["0"]]),decreasing=T) ]

 PKR.diff <-
  diffMat(
   As=paste(c("emb_Prkra_Mut_r1","emb_Prkra_Mut_r2","emb_Prkra_Mut_r3"),".CPs.normCent.mat",sep=""),
   Bs=paste(c("emb_Prkra_WT_r1", "emb_Prkra_WT_r1"),                    ".CPs.normCent.mat",sep=""),
   FCE="median"
  )
# >> order from -1 => +1 at "0" position ( == increasing )
 PKR.diff   <- PKR.diff[   order(   (PKR.diff[,"0"]),   decreasing=F) , ]
 tmp.dt     <- data.table( PKR.diff,keep.rownames="ID" )
 PKR.CPs.dt <- merge( PKR.CPs.dt,tmp.dt,by="ID",all=T )
# >> order from biggest changes (abs.values) to smallest at "0" position ( == decreasing )
 PKR.CPs.dt <- PKR.CPs.dt[ order(abs(PKR.CPs.dt[["0"]]),decreasing=T) ]

 TAR.diff <-
  diffMat(
   As=paste(c("emb_Tarbp2_Mut_r1","emb_Tarbp2_Mut_r2","emb_Tarbp2_Mut_r3"),".CPs.normCent.mat",sep=""),
   Bs=paste(c("emb_Tarbp2_WT_r1", "emb_Tarbp2_WT_r2", "emb_Tarbp2_WT_r3"), ".CPs.normCent.mat",sep=""),
   FCE="median"
  )
# >> order from -1 => +1 at "0" position ( == increasing )
 TAR.diff   <- TAR.diff[   order(   (TAR.diff[,"0"]),   decreasing=F) , ]
 tmp.dt     <- data.table( TAR.diff,keep.rownames="ID" )
 TAR.CPs.dt <- merge( TAR.CPs.dt,tmp.dt,by="ID",all=T )
# >> order from biggest changes (abs.values) to smallest at "0" position ( == decreasing )
 TAR.CPs.dt <- TAR.CPs.dt[ order(abs(TAR.CPs.dt[["0"]]),decreasing=T) ]

 GNT.diff <-
  diffMat(
   As=paste(c("x3_GNT_homozygot","x4_GNT_homozygot","x9_GNT_homozygot"),".CPs.normCent.mat",sep=""),
   Bs=paste(c("x11_WT",          "x14_WT",          "x16_WT"),          ".CPs.normCent.mat",sep=""),
   FCE="median"
  )
# >> order from -1 => +1 at "0" position ( == increasing )
 GNT.diff   <- GNT.diff[   order(   (GNT.diff[,"0"]),   decreasing=F) , ]
 tmp.dt     <- data.table( GNT.diff,keep.rownames="ID" )
 GNT.CPs.dt <- merge( GNT.CPs.dt,tmp.dt,by="ID",all=T )
# >> order from biggest changes (abs.values) to smallest at "0" position ( == decreasing )
 GNT.CPs.dt <- GNT.CPs.dt[ order(abs(GNT.CPs.dt[["0"]]),decreasing=T) ]

if ( T ) {

# (1) EMB-ES
 rLAB <- EMB.CPs.dt[ , paste0(" ",Name," ") ]; names(rLAB) <- EMB.CPs.dt[,ID]
 pheatmap(
  EMB.diff[ rownames(EMB.diff) %in% head( EMB.CPs.dt[ det.STR=="5p" , ID ],TOPN ) , ],cluster_rows=F,cluster_cols=F,
  main=paste("EMB-miRNA-5p : top ",TOPN,sep=""),border_color="black",cellwidth=50,cellheight=50,
  color=colorRampPalette(c("blue","white","red"))( length(breakList) ),breaks=breakList,
  labels_row=rLAB[rownames(head(EMB.diff[ rownames(EMB.diff) %in% EMB.CPs.dt[ det.STR=="5p" , ID ] , ],TOPN))],
  filename=paste("EMB-fidelity-5p-top",TOPN,"-",ONAME,"-",pepDate(),".pdf",sep=""),width=25,height=040,fontsize=50
 )
 saveRDS( EMB.diff[ rownames(EMB.diff) %in% head( EMB.CPs.dt[ det.STR=="5p" , ID ],TOPN ) , ],file=paste0("EMB-5p.top",TOPN,"-",pepDate(),".rda") )
 pheatmap(
  EMB.diff[ rownames(EMB.diff) %in% head( EMB.CPs.dt[ det.STR=="3p" , ID ],TOPN ) , ],cluster_rows=F,cluster_cols=F,
  main=paste("EMB-miRNA-3p : top ",TOPN,sep=""),border_color="black",cellwidth=50,cellheight=50,
  color=colorRampPalette(c("blue","white","red"))( length(breakList) ),breaks=breakList,
  labels_row=rLAB[rownames(head(EMB.diff[ rownames(EMB.diff) %in% EMB.CPs.dt[ det.STR=="3p" , ID ] , ],TOPN))],
  filename=paste("EMB-fidelity-3p-top",TOPN,"-",ONAME,"-",pepDate(),".pdf",sep=""),width=25,height=040,fontsize=50
 )
 saveRDS( EMB.diff[ rownames(EMB.diff) %in% head( EMB.CPs.dt[ det.STR=="3p" , ID ],TOPN ) , ],file=paste0("EMB-3p.top",TOPN,"-",pepDate(),".rda") )

# (2) ESC-ES
 rLAB <- ESC.CPs.dt[ , paste0(" ",Name," ") ]; names(rLAB) <- ESC.CPs.dt[,ID]
 pheatmap(
  ESC.diff[ rownames(ESC.diff) %in% head( ESC.CPs.dt[ det.STR=="5p" , ID ],TOPN ) , ],cluster_rows=F,cluster_cols=F,
  main=paste("ESC-miRNA-5p : top ",TOPN,sep=""),border_color="black",cellwidth=50,cellheight=50,
  color=colorRampPalette(c("blue","white","red"))( length(breakList) ),breaks=breakList,
  labels_row=rLAB[rownames(head(ESC.diff[ rownames(ESC.diff) %in% ESC.CPs.dt[ det.STR=="5p" , ID ] , ],TOPN))],
  filename=paste("ESC-fidelity-5p-top",TOPN,"-",ONAME,"-",pepDate(),".pdf",sep=""),width=25,height=040,fontsize=50
 )
 saveRDS( ESC.diff[ rownames(ESC.diff) %in% head( ESC.CPs.dt[ det.STR=="5p" , ID ],TOPN ) , ],file=paste0("ESC-5p.top",TOPN,"-",pepDate(),".rda") )
 pheatmap(
  ESC.diff[ rownames(ESC.diff) %in% head( ESC.CPs.dt[ det.STR=="3p" , ID ],TOPN ) , ],cluster_rows=F,cluster_cols=F,
  main=paste("ESC-miRNA-3p : top ",TOPN,sep=""),border_color="black",cellwidth=50,cellheight=50,
  color=colorRampPalette(c("blue","white","red"))( length(breakList) ),breaks=breakList,
  labels_row=rLAB[rownames(head(ESC.diff[ rownames(ESC.diff) %in% ESC.CPs.dt[ det.STR=="3p" , ID ] , ],TOPN))],
  filename=paste("ESC-fidelity-3p-top",TOPN,"-",ONAME,"-",pepDate(),".pdf",sep=""),width=25,height=040,fontsize=50
 )
 saveRDS( ESC.diff[ rownames(ESC.diff) %in% head( ESC.CPs.dt[ det.STR=="3p" , ID ],TOPN ) , ],file=paste0("ESC-3p.top",TOPN,"-",pepDate(),".rda") )
 
# (3) PKR-Pu
 rLAB <- PKR.CPs.dt[ , paste0(" ",Name," ") ]; names(rLAB) <- PKR.CPs.dt[,ID]
 pheatmap(
  PKR.diff[ rownames(PKR.diff) %in% head( PKR.CPs.dt[ det.STR=="5p" , ID ],TOPN ) , ],cluster_rows=F,cluster_cols=F,
  main=paste("PKR-miRNA-5p : top ",TOPN,sep=""),border_color="black",cellwidth=50,cellheight=50,
  color=colorRampPalette(c("blue","white","red"))( length(breakList) ),breaks=breakList,
  labels_row=rLAB[rownames(head(PKR.diff[ rownames(PKR.diff) %in% PKR.CPs.dt[ det.STR=="5p" , ID ] , ],TOPN))],
  filename=paste("PKR-fidelity-5p-top",TOPN,"-",ONAME,"-",pepDate(),".pdf",sep=""),width=25,height=040,fontsize=50
 )
 saveRDS( PKR.diff[ rownames(PKR.diff) %in% head( PKR.CPs.dt[ det.STR=="5p" , ID ],TOPN ) , ],file=paste0("PKR-5p.top",TOPN,"-",pepDate(),".rda") )
 pheatmap(
  PKR.diff[ rownames(PKR.diff) %in% head( PKR.CPs.dt[ det.STR=="3p" , ID ],TOPN ) , ],cluster_rows=F,cluster_cols=F,
  main=paste("PKR-miRNA-3p : top ",TOPN,sep=""),border_color="black",cellwidth=50,cellheight=50,
  color=colorRampPalette(c("blue","white","red"))( length(breakList) ),breaks=breakList,
  labels_row=rLAB[rownames(head(PKR.diff[ rownames(PKR.diff) %in% PKR.CPs.dt[ det.STR=="3p" , ID ] , ],TOPN))],
  filename=paste("PKR-fidelity-3p-top",TOPN,"-",ONAME,"-",pepDate(),".pdf",sep=""),width=25,height=040,fontsize=50
 )
 saveRDS( PKR.diff[ rownames(PKR.diff) %in% head( PKR.CPs.dt[ det.STR=="3p" , ID ],TOPN ) , ],file=paste0("PKR-3p.top",TOPN,"-",pepDate(),".rda") )
 
# (4) TAR-Pu
 rLAB <- TAR.CPs.dt[ , paste0(" ",Name," ") ]; names(rLAB) <- TAR.CPs.dt[,ID]
 pheatmap(
  TAR.diff[ rownames(TAR.diff) %in% head( TAR.CPs.dt[ det.STR=="5p" , ID ],TOPN ) , ],cluster_rows=F,cluster_cols=F,
  main=paste("TAR-miRNA-5p : top ",TOPN,sep=""),border_color="black",cellwidth=50,cellheight=50,
  color=colorRampPalette(c("blue","white","red"))( length(breakList) ),breaks=breakList,
  labels_row=rLAB[rownames(head(TAR.diff[ rownames(TAR.diff) %in% TAR.CPs.dt[ det.STR=="5p" , ID ] , ],TOPN))],
  filename=paste("TAR-fidelity-5p-top",TOPN,"-",ONAME,"-",pepDate(),".pdf",sep=""),width=25,height=040,fontsize=50
 )
 saveRDS( TAR.diff[ rownames(TAR.diff) %in% head( TAR.CPs.dt[ det.STR=="5p" , ID ],TOPN ) , ],file=paste0("TAR-5p.top",TOPN,"-",pepDate(),".rda") )
 pheatmap(
  TAR.diff[ rownames(TAR.diff) %in% head( TAR.CPs.dt[ det.STR=="3p" , ID ],TOPN ) , ],cluster_rows=F,cluster_cols=F,
  main=paste("TAR-miRNA-3p : top ",TOPN,sep=""),border_color="black",cellwidth=50,cellheight=50,
  color=colorRampPalette(c("blue","white","red"))( length(breakList) ),breaks=breakList,
  labels_row=rLAB[rownames(head(TAR.diff[ rownames(TAR.diff) %in% TAR.CPs.dt[ det.STR=="3p" , ID ] , ],TOPN))],
  filename=paste("TAR-fidelity-3p-top",TOPN,"-",ONAME,"-",pepDate(),".pdf",sep=""),width=25,height=040,fontsize=50
 )
 saveRDS( TAR.diff[ rownames(TAR.diff) %in% head( TAR.CPs.dt[ det.STR=="3p" , ID ],TOPN ) , ],file=paste0("TAR-3p.top",TOPN,"-",pepDate(),".rda") )
 
# (5) GNT-VB
 rLAB <- GNT.CPs.dt[ , paste0(" ",Name," ") ]; names(rLAB) <- GNT.CPs.dt[,ID]
 pheatmap(
  GNT.diff[ rownames(GNT.diff) %in% head( GNT.CPs.dt[ det.STR=="5p" , ID ],TOPN ) , ],cluster_rows=F,cluster_cols=F,
  main=paste("GNT-miRNA-5p : top ",TOPN,sep=""),border_color="black",cellwidth=50,cellheight=50,
  color=colorRampPalette(c("blue","white","red"))( length(breakList) ),breaks=breakList,
  labels_row=rLAB[rownames(head(GNT.diff[ rownames(GNT.diff) %in% GNT.CPs.dt[ det.STR=="5p" , ID ] , ],TOPN))],
  filename=paste("GNT-fidelity-5p-top",TOPN,"-",ONAME,"-",pepDate(),".pdf",sep=""),width=25,height=040,fontsize=50
 )
 saveRDS( GNT.diff[ rownames(GNT.diff) %in% head( GNT.CPs.dt[ det.STR=="5p" , ID ],TOPN ) , ],file=paste0("GNT-5p.top",TOPN,"-",pepDate(),".rda") )
 pheatmap(
  GNT.diff[ rownames(GNT.diff) %in% head( GNT.CPs.dt[ det.STR=="3p" , ID ],TOPN ) , ],cluster_rows=F,cluster_cols=F,
  main=paste("GNT-miRNA-3p : top ",TOPN,sep=""),border_color="black",cellwidth=50,cellheight=50,
  color=colorRampPalette(c("blue","white","red"))( length(breakList) ),breaks=breakList,
  labels_row=rLAB[rownames(head(GNT.diff[ rownames(GNT.diff) %in% GNT.CPs.dt[ det.STR=="3p" , ID ] , ],TOPN))],
  filename=paste("GNT-fidelity-3p-top",TOPN,"-",ONAME,"-",pepDate(),".pdf",sep=""),width=25,height=040,fontsize=50
 )
 saveRDS( GNT.diff[ rownames(GNT.diff) %in% head( GNT.CPs.dt[ det.STR=="3p" , ID ],TOPN ) , ],file=paste0("GNT-3p.top",TOPN,"-",pepDate(),".rda") )
 
}

if ( F ) {
EMB.CPs.dt[["DicerX_WT_6_r1"]]  <- NULL
EMB.CPs.dt[["DicerX_WT_8B_r2"]] <- NULL
EMB.CPs.dt[["DicerX_WT_1_r3"]]  <- NULL
EMB.CPs.dt[["CP"]]                      <- NULL
ESC.CPs.dt[["RS7_Mos_r1"]]            <- NULL
ESC.CPs.dt[["RS7_Mos_r2"]]            <- NULL
ESC.CPs.dt[["RS7_Mos_r3"]]            <- NULL
ESC.CPs.dt[["CP"]]                      <- NULL
xwb     <- createWorkbook()
addWorksheet(   wb=xwb,sheetName="EMB" )
writeDataTable( wb=xwb,sheet=    "EMB",x=EMB.CPs.dt,colNames=T,rowNames=F )
addWorksheet(   wb=xwb,sheetName="ESC" )
writeDataTable( wb=xwb,sheet=    "ESC",x=ESC.CPs.dt,colNames=T,rowNames=F )
saveWorkbook(   wb=xwb,file=paste("fidelity-",pepDate(),".xlsx",sep=""),overwrite=T )
}

}

