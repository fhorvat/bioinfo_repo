library(data.table)
library(GenomicAlignments)
library(openxlsx)

load( "mirAnnot.dt.rda",verbose=T )

#= FUN
cleavageEfficiency <- function( CRDS,PRELOC,CRDDEV=2,xGR=tmp.gr,PLOT=T,TIT=NULL ) {

 if ( length(CRDS)!=4 ) {
  cat( "\n ! length(CRDS) = ",length(CRDS)," ! \n ! \"CRDS\" musts contain 4 values : \"Drosha-5p-CP\", \"Dicer-5p-CP\", \"Dicer-3p-CP\", \"Drosha-3p-CP\" ! \n\n",sep="" )
  stop()
 }

 if ( is.null(TIT) ) {
  TIT=PRELOC
 }

 if ( sub( "^.*[:]","",as.character(PRELOC) )=="+" ) {
  CRDS   <- sort( x=as.numeric(CRDS),decreasing=F )
 } else                                              {
  CRDS   <- sort( x=as.numeric(CRDS),decreasing=T )
 }

 LOCS <-
 list(
  mi5p = apply( X=as.matrix(x=CJ( seq(CRDS[1]-CRDDEV,CRDS[1]+CRDDEV),seq(CRDS[2]-CRDDEV,CRDS[2]+CRDDEV),unique=T )),MARGIN=1,FUN=paste,collapse="-" ),
  loop = apply( X=as.matrix(x=CJ( seq(CRDS[2]-CRDDEV,CRDS[2]+CRDDEV),seq(CRDS[3]-CRDDEV,CRDS[3]+CRDDEV),unique=T )),MARGIN=1,FUN=paste,collapse="-" ),
  mi3p = apply( X=as.matrix(x=CJ( seq(CRDS[3]-CRDDEV,CRDS[3]+CRDDEV),seq(CRDS[4]-CRDDEV,CRDS[4]+CRDDEV),unique=T )),MARGIN=1,FUN=paste,collapse="-" ),
  ov5p = apply( X=as.matrix(x=CJ( seq(CRDS[1]-CRDDEV,CRDS[1]+CRDDEV),seq(CRDS[3]-CRDDEV,CRDS[3]+CRDDEV),unique=T )),MARGIN=1,FUN=paste,collapse="-" ),
  ov3p = apply( X=as.matrix(x=CJ( seq(CRDS[2]-CRDDEV,CRDS[2]+CRDDEV),seq(CRDS[4]-CRDDEV,CRDS[4]+CRDDEV),unique=T )),MARGIN=1,FUN=paste,collapse="-" ),
  preL = apply( X=as.matrix(x=CJ( seq(CRDS[1]-CRDDEV,CRDS[1]+CRDDEV),seq(CRDS[4]-CRDDEV,CRDS[4]+CRDDEV),unique=T )),MARGIN=1,FUN=paste,collapse="-" )
 )

 pre.gr <- GenomicRanges::findOverlaps( query=xGR,subject=GRanges(PRELOC),type="within",ignore.strand=F )
 pre.gr <- xGR[ unique(queryHits(pre.gr)) ]

 if ( sub( "^.*[:]","",as.character(PRELOC) )=="+" ) {
  cnt.dt <- data.table( reads=paste0( start(pre.gr),"-",  end(pre.gr) ) )[ , .N , by="reads" ]
 } else                                              {
  cnt.dt <- data.table( reads=paste0(   end(pre.gr),"-",start(pre.gr) ) )[ , .N , by="reads" ]
 }

 TOTN   <- cnt.dt[ , sum(N) ]

   perc.mi5p = ( sum( cnt.dt[ reads %in% LOCS[["mi5p"]] ][["N"]] ) / TOTN )*100
   perc.loop = ( sum( cnt.dt[ reads %in% LOCS[["loop"]] ][["N"]] ) / TOTN )*100
   perc.mi3p = ( sum( cnt.dt[ reads %in% LOCS[["mi3p"]] ][["N"]] ) / TOTN )*100
   perc.ov5p = ( sum( cnt.dt[ reads %in% LOCS[["ov5p"]] ][["N"]] ) / TOTN )*100
   perc.ov3p = ( sum( cnt.dt[ reads %in% LOCS[["ov3p"]] ][["N"]] ) / TOTN )*100
   perc.preL = ( sum( cnt.dt[ reads %in% LOCS[["preL"]] ][["N"]] ) / TOTN )*100

 out.dt <-
  data.table(
   mi5p = perc.mi5p,
   loop = perc.loop,
   mi3p = perc.mi3p,
   ov5p = perc.ov5p,
   ov3p = perc.ov3p,
   preL = perc.preL,
   totN = TOTN
  )

 XLIM=c(1,70)
 XCRD=c(00,10,30,40,60,70)
 YLIM=c(-1,7)
 XL <- c( XCRD[2],XCRD[3],XCRD[4],XCRD[2],XCRD[3],XCRD[2],XCRD[1],XCRD[2],XCRD[3],XCRD[4],XCRD[5] )
 XR <- c( XCRD[3],XCRD[4],XCRD[5],XCRD[4],XCRD[5],XCRD[5],XCRD[2],XCRD[3],XCRD[4],XCRD[5],XCRD[6] )
 YB <- rep.int( c(5.05,4.05,3.05,2.05,1.25,1.15,1.25,1.15,1.25),c(3,1,1,1,1,1,1,1,1) )
 YT <- YB + rep.int( c(0.90,0.50,0.70,0.50,0.70,0.50),c(6,1,1,1,1,1) )
 if ( PLOT ) {

  par( bg="white",xpd=T )
  plot( -1000,-1000,xlab="",ylab="",xaxt="n",yaxt="n",main=TIT,xlim=XLIM,ylim=YLIM,bty="n" )
  rect(
   xleft   = XL,
   xright  = XR,
   ybottom = YB,
   ytop    = YT,
   border  = "black",
   col     = rep.int( c("gray50","gray50","blue","black","red","gray50"),c(6,1,1,1,1,1) )
  )
  text(
   x      = ((XL+XR)/2),
   y      = ((YB+YT)/2) - rep.int( c(0,1,1,1,1,1),c(6,1,1,1,1,1) ),
   labels = c( sprintf( "%4.1f%%",unlist(out.dt)[ c("mi5p","loop","mi3p","ov5p","ov3p","preL") ] ),c("","5p","loop","3p","") ),
   font   = 2
  )
  points(
   x=XCRD[3:4],y=c(YB[length(YB)]-00.5,YB[length(YB)]-00.5),pch=24,bg="cyan"
  )
  text(
   x=XCRD[3:4],y=c(YB[length(YB)]-01.5,YB[length(YB)]-01.5),labels=c( "3' 5p end\ncleavage\npoint","5' 3p end\ncleavage\npoint" ),font=2
  )
  text(
   x      = sum(XCRD[4:5])/2,
   y      = YLIM[2]*0.95,
   labels = sprintf( "Sum of visualized reads = %4.1f%%",sum(unlist(out.dt)[ c("mi5p","loop","mi3p","ov5p","ov3p","preL") ]) ),
   font   = 2
  )
  text(
   x      = sum(XCRD[4:5])/2,
   y      = YLIM[2]*0.90,
#   labels = sprintf( "Total number of reads = %s",numSepPointsFormat( n=out.dt[["totN"]] ) ),
   labels = sprintf( "Total number of reads = %s",out.dt[["totN"]] ),
   font   = 2
  )

 }

 out.dt[["name"]] <- TIT
 out.dt[["gLoc"]] <- PRELOC
 out.dt           <- out.dt[ , c("name","gLoc","mi5p","loop","mi3p","ov5p","ov3p","preL","totN") , with=F ]

 return( out.dt )

}

plot_cEff <- function( par.dt=cEff.dt,WT,KO,MFROW=NULL,TIT=NULL,YLABS=setNames( c("wt","DicerX"),c("WT","KO") ) ) {

 if ( length(unique(par.dt[["name"]]))!=1 ) {
  stop("\n\n !!! Multiple miRNAs in the input \"par.dt\" !!!\n\n")
 }

 if ( is.null(MFROW) ) {
  MFROW <- c( 2,max(c(length(WT),length(KO))) )
 }

 par.dt <- par.dt[ sample %in% c(WT,KO) ]

 XLIM=c(1,70)
 XCRD=c(00,10,30,40,60,70)
 YLIM=c(-1,7)
 XL <- c( XCRD[2],XCRD[3],XCRD[4],XCRD[2],XCRD[3],XCRD[2],XCRD[1],XCRD[2],XCRD[3],XCRD[4],XCRD[5] )
 XR <- c( XCRD[3],XCRD[4],XCRD[5],XCRD[4],XCRD[5],XCRD[5],XCRD[2],XCRD[3],XCRD[4],XCRD[5],XCRD[6] )
 YB <- rep.int( c(5.05,4.05,3.05,2.05,1.25,1.15,1.25,1.15,1.25),c(3,1,1,1,1,1,1,1,1) )
 YT <- YB + rep.int( c(0.90,0.50,0.70,0.50,0.70,0.50),c(6,1,1,1,1,1) )

 par( bg="white",xpd=T,mfrow=MFROW )
 for ( CON in c("WT","KO") ) {
 for ( i in seq(max(c(length(WT),length(KO)))) ) {
 if ( is.na( get(CON)[i] ) ) {
  frame()
 } else          {
 if ( is.null(TIT) ) {
  TIT <- paste0( get(CON)[i]," : ",unique(par.dt[["name"]]) )
 }
 cat("   > ",get(CON)[i],"\n",sep="")
 plot( -1000,-1000,xlab="",ylab="",xaxt="n",yaxt="n",main=TIT,xlim=XLIM,ylim=YLIM,bty="n" )
# if ( i==1 ) { title( ylab=YLABS[CON],font=2 ) }
 if ( i==1 ) { text( x=-3.00,y=mean(c(YB[length(YB)-1],YT[1])),labels=YLABS[CON],font=2,srt=90,cex=2.00 ) }
 rect(
  xleft   = XL,
  xright  = XR,
  ybottom = YB,
  ytop    = YT,
  border  = "black",
  col     = rep.int( c("gray50","gray50","blue","black","red","gray50"),c(6,1,1,1,1,1) )
 )
 text(
  x      = ((XL+XR)/2),
  y      = ((YB+YT)/2) - rep.int( c(0,1,1,1,1,1),c(6,1,1,1,1,1) ),
  labels = c( sprintf( "%6.3f%%",unlist(par.dt[ sample==get(CON)[i],c("mi5p","loop","mi3p","ov5p","ov3p","preL"),with=F ]) ),c("","5p","loop","3p","") ),
  font   = 2
 )
 points(
  x=XCRD[3:4],y=c(YB[length(YB)]-00.5,YB[length(YB)]-00.5),pch=24,bg="cyan"
 )
 text(
  x=XCRD[3:4],y=c(YB[length(YB)]-01.5,YB[length(YB)]-01.5),labels=c( "3' 5p end\ncleavage\npoint","5' 3p end\ncleavage\npoint" ),font=2
 )
 text(
  x      = sum(XCRD[4:5])/2,
  y      = YLIM[2]*0.95,
  labels = sprintf( "Sum of visualized reads = %6.3f%%",sum(unlist(par.dt[ sample==get(CON)[i],c("mi5p","loop","mi3p","ov5p","ov3p","preL"),with=F ])) ),
  font   = 2
 )
 text(
  x      = sum(XCRD[4:5])/2,
  y      = YLIM[2]*0.90,
#  labels = sprintf( "Total number of reads = %s",numSepPointsFormat( n=par.dt[ sample==get(CON)[i] ][["totN"]] ) ),
  labels = sprintf( "Total number of reads = %s",par.dt[ sample==get(CON)[i] ][["totN"]] ),
  font   = 2
 )
 TIT <- NULL
 }
 }
 }

 return( TRUE )

}

avrPlot_cEff <- function( par.dt=cEff.dt,WT,KO,MFROW=c(1,2),YLABS=setNames( c("wt","DicerX"),c("WT","KO") ) ) {

 if ( length(unique(par.dt[["name"]]))!=1 ) {
  stop("\n\n !!! Multiple miRNAs in the input \"par.dt\" !!!\n\n")
 }

 par.dt <- par.dt[ sample %in% c(WT,KO) ]

 avr.dt <- data.table( name=rep.int(unique(par.dt[["name"]]),2),gLoc=rep.int(unique(par.dt[["gLoc"]]),2),sample=as.vector(YLABS) )
 for ( iCLN in c("mi5p","loop","mi3p","ov5p","ov3p","preL","totN") ) {
  avr.dt[[iCLN]]                                <- as.numeric("")
  avr.dt[ sample==as.vector(YLABS[1]) ][[iCLN]] <- mean( par.dt[ sample %in% WT ][[iCLN]] )
  avr.dt[ sample==as.vector(YLABS[2]) ][[iCLN]] <- mean( par.dt[ sample %in% KO ][[iCLN]] )
 }
 par.dt <- avr.dt

 TIT=NULL
 XLIM=c(1,70)
 XCRD=c(00,10,30,40,60,70)
 YLIM=c(-1,7)
 XL <- c( XCRD[2],XCRD[3],XCRD[4],XCRD[2],XCRD[3],XCRD[2],XCRD[1],XCRD[2],XCRD[3],XCRD[4],XCRD[5] )
 XR <- c( XCRD[3],XCRD[4],XCRD[5],XCRD[4],XCRD[5],XCRD[5],XCRD[2],XCRD[3],XCRD[4],XCRD[5],XCRD[6] )
 YB <- rep.int( c(5.05,4.05,3.05,2.05,1.25,1.15,1.25,1.15,1.25),c(3,1,1,1,1,1,1,1,1) )
 YT <- YB + rep.int( c(0.90,0.50,0.70,0.50,0.70,0.50),c(6,1,1,1,1,1) )

 par( bg="white",xpd=T,mfrow=MFROW )
 for ( CON in c("WT","KO") ) {
 if ( is.null(TIT) ) {
  TIT <- paste0( CON," : ",unique(par.dt[["name"]]) )
 }
 cat("   > ",CON,"\n",sep="")
 plot( -1000,-1000,xlab="",ylab="",xaxt="n",yaxt="n",main=TIT,xlim=XLIM,ylim=YLIM,bty="n" )
# text( x=-3.00,y=mean(c(YB[length(YB)-1],YT[1])),labels=YLABS[CON],font=2,srt=90,cex=2.00 )
 rect(
  xleft   = XL,
  xright  = XR,
  ybottom = YB,
  ytop    = YT,
  border  = "black",
  col     = rep.int( c("gray50","gray50","blue","black","red","gray50"),c(6,1,1,1,1,1) )
 )
 text(
  x      = ((XL+XR)/2),
  y      = ((YB+YT)/2) - rep.int( c(0,1,1,1,1,1),c(6,1,1,1,1,1) ),
  labels = c( sprintf( "%6.3f%%",unlist(par.dt[ sample==YLABS[CON],c("mi5p","loop","mi3p","ov5p","ov3p","preL"),with=F ]) ),c("","5p","loop","3p","") ),
  font   = 2
 )
 points(
  x=XCRD[3:4],y=c(YB[length(YB)]-00.5,YB[length(YB)]-00.5),pch=24,bg="cyan"
 )
 text(
  x=XCRD[3:4],y=c(YB[length(YB)]-01.5,YB[length(YB)]-01.5),labels=c( "3' 5p end\ncleavage\npoint","5' 3p end\ncleavage\npoint" ),font=2
 )
 text(
  x      = sum(XCRD[4:5])/2,
  y      = YLIM[2]*0.95,
  labels = sprintf( "Sum of visualized reads = %6.3f%%",sum(unlist(par.dt[ sample==YLABS[CON],c("mi5p","loop","mi3p","ov5p","ov3p","preL"),with=F ])) ),
  font   = 2
 )
 text(
  x      = sum(XCRD[4:5])/2,
  y      = YLIM[2]*0.90,
#  labels = sprintf( "Total number of reads = %s",numSepPointsFormat( n=round(par.dt[ sample==YLABS[CON] ][["totN"]]) ) ),
  labels = sprintf( "Total number of reads = %s",round(par.dt[ sample==YLABS[CON] ][["totN"]],1) ),
  font   = 2
 )
 TIT <- NULL
 }

 return( avr.dt )

}

#=

#= miR-15a  : chr14:61632038-61632096:- | c( 61632075,61632096,61632038,61632059 )
#= miR-7068 : chr8:72470016-72470089:-  | c( 72470069,72470089,72470016,72470036 )
#= miR-145a : chr18:61647827-61647888:- | c( 61647866,61647888,61647827,61647848 )
inp.mirna.dt <-
 data.table(
  N = c("miR-15a","miR-7068","miR-145a"),
  L = c("chr14:61632038-61632096:-","chr8:72470016-72470089:-","chr18:61647827-61647888:-"),
  C = c("61632075,61632096,61632038,61632059","72470069,72470089,72470016,72470036","61647866,61647888,61647827,61647848")
 )

if (F)   {

cEff.dt <- data.table()

#= (1) EMB-15.5
STOR="/storage/brno1-cerit/home/pepap/Taborska_smallRNAlibs_DicerXXE_191217/BAM/"
SUFF=".se.Aligned.sortedByCoord.out.bam"
BAMFILES <- system( command=paste0("/bin/ls ",STOR,"*/*",SUFF),intern=T )
OUTNAMES <- gsub( "[0-9][0-9][.]","",gsub( "[/].*$","",gsub( STOR,"",system( command=paste0("/bin/ls ",STOR,"*/*",SUFF),intern=T ) ) ) )

cat( "\n",sep="" )
for ( ibf in seq_along(BAMFILES) ) {

 cat( " => ",OUTNAMES[ibf],"\n",sep="" )

 tmp.tga <- readGAlignments( file=BAMFILES[ibf],param=ScanBamParam( tag=c("HI","NH","nM") ) )
 tmp.tga <- tmp.tga[ njunc(tmp.tga)==0 ]
 tmp.gr  <- as( object=tmp.tga,Class="GRanges" )

 for ( imir in seq_along(inp.mirna.dt[["N"]]) ) {
  cat( "  >> ",inp.mirna.dt[["N"]][ imir ],"\n",sep="" )
  dir.create( path=paste0("PDF/",inp.mirna.dt[["N"]][ imir ]),showWarnings=F,recursive=T )
  pdf( file=paste0("PDF/",inp.mirna.dt[["N"]][ imir ],"/",OUTNAMES[ibf],"-",pepDate(),".pdf"),width=25,height=10,pointsize=15 )
  tmp.dt  <-
   cleavageEfficiency(
    CRDS   = unlist(strsplit( x=inp.mirna.dt[["C"]][ imir ],split=",",fixed=T )),
    PRELOC = inp.mirna.dt[["L"]][ imir ],
    CRDDEV = 2,
    xGR    = tmp.gr,
    PLOT   = T,
    TIT    = inp.mirna.dt[["N"]][ imir ]
   )
  dev.off()
  tmp.dt[["sample"]] <- OUTNAMES[ibf]
  cEff.dt <- rbind( cEff.dt,tmp.dt )
 }

}

#= (2) ESC
STOR="/storage/brno1-cerit/home/pepap/Eliska.Taborska/BAM/"
SUFF=".se.Aligned.sortedByCoord.out.bam"
BAMFILES <- system( command=paste0("/bin/ls ",STOR,"*/*",SUFF),intern=T )
OUTNAMES <- gsub( "[0-9][0-9][.]","",gsub( "[/].*$","",gsub( STOR,"",system( command=paste0("/bin/ls ",STOR,"*/*",SUFF),intern=T ) ) ) )

cat( "\n",sep="" )
for ( ibf in seq_along(BAMFILES) ) {

 cat( " => ",OUTNAMES[ibf],"\n",sep="" )

 tmp.tga <- readGAlignments( file=BAMFILES[ibf],param=ScanBamParam( tag=c("HI","NH","nM") ) )
 tmp.tga <- tmp.tga[ njunc(tmp.tga)==0 ]
 tmp.gr  <- as( object=tmp.tga,Class="GRanges" )

 for ( imir in seq_along(inp.mirna.dt[["N"]]) ) {
  cat( "  >> ",inp.mirna.dt[["N"]][ imir ],"\n",sep="" )
  dir.create( path=paste0("PDF/",inp.mirna.dt[["N"]][ imir ]),showWarnings=F,recursive=T )
  pdf( file=paste0("PDF/",inp.mirna.dt[["N"]][ imir ],"/",OUTNAMES[ibf],"-",pepDate(),".pdf"),width=25,height=10,pointsize=15 )
  tmp.dt  <-
   cleavageEfficiency(
    CRDS   = unlist(strsplit( x=inp.mirna.dt[["C"]][ imir ],split=",",fixed=T )),
    PRELOC = inp.mirna.dt[["L"]][ imir ],
    CRDDEV = 2,
    xGR    = tmp.gr,
    PLOT   = T,
    TIT    = inp.mirna.dt[["N"]][ imir ]
   )
  dev.off()
  tmp.dt[["sample"]] <- OUTNAMES[ibf]
  cEff.dt <- rbind( cEff.dt,tmp.dt )
 }

}

#= (3) Pullagura : PKR + TAR
STOR="/storage/brno1-cerit/home/pepap/Pullagura_2018_Genetics_PRJNA423238/BAM/READs/"
SUFF=".se.Aligned.sortedByCoord.out.bam"
BAMFILES <- system( command=paste0("/bin/ls ",STOR,"*/*",SUFF),intern=T )
OUTNAMES <- gsub( "[0-9][0-9][.]","",gsub( "[/].*$","",gsub( STOR,"",system( command=paste0("/bin/ls ",STOR,"*/*",SUFF),intern=T ) ) ) )

cat( "\n",sep="" )
for ( ibf in seq_along(BAMFILES) ) {

 cat( " => ",OUTNAMES[ibf],"\n",sep="" )

 tmp.tga <- readGAlignments( file=BAMFILES[ibf],param=ScanBamParam( tag=c("HI","NH","nM") ) )
 tmp.tga <- tmp.tga[ njunc(tmp.tga)==0 ]
 tmp.gr  <- as( object=tmp.tga,Class="GRanges" )

 for ( imir in seq_along(inp.mirna.dt[["N"]]) ) {
  cat( "  >> ",inp.mirna.dt[["N"]][ imir ],"\n",sep="" )
  dir.create( path=paste0("PDF/",inp.mirna.dt[["N"]][ imir ]),showWarnings=F,recursive=T )
  pdf( file=paste0("PDF/",inp.mirna.dt[["N"]][ imir ],"/",OUTNAMES[ibf],"-",pepDate(),".pdf"),width=25,height=10,pointsize=15 )
  tmp.dt  <-
   cleavageEfficiency(
    CRDS   = unlist(strsplit( x=inp.mirna.dt[["C"]][ imir ],split=",",fixed=T )),
    PRELOC = inp.mirna.dt[["L"]][ imir ],
    CRDDEV = 2,
    xGR    = tmp.gr,
    PLOT   = T,
    TIT    = inp.mirna.dt[["N"]][ imir ]
   )
  dev.off()
  tmp.dt[["sample"]] <- OUTNAMES[ibf]
  cEff.dt <- rbind( cEff.dt,tmp.dt )
 }

}

#= (4) GNT
STOR="/storage/brno1-cerit/home/pepap/Valeria.Buccheri/01.Dicer-helicase-domain-mutations/BAM.MERGED/"
SUFF=".se.Aligned.sortedByCoord.out.bam"
BAMFILES <- system( command=paste0("/bin/ls ",STOR,"*/*",SUFF),intern=T )
OUTNAMES <- gsub( "[0-9][0-9][.]","",gsub( "[/].*$","",gsub( STOR,"",system( command=paste0("/bin/ls ",STOR,"*/*",SUFF),intern=T ) ) ) )

cat( "\n",sep="" )
for ( ibf in seq_along(BAMFILES) ) {

 cat( " => ",OUTNAMES[ibf],"\n",sep="" )

 tmp.tga <- readGAlignments( file=BAMFILES[ibf],param=ScanBamParam( tag=c("HI","NH","nM") ) )
 tmp.tga <- tmp.tga[ njunc(tmp.tga)==0 ]
 tmp.gr  <- as( object=tmp.tga,Class="GRanges" )

 for ( imir in seq_along(inp.mirna.dt[["N"]]) ) {
  cat( "  >> ",inp.mirna.dt[["N"]][ imir ],"\n",sep="" )
  dir.create( path=paste0("PDF/",inp.mirna.dt[["N"]][ imir ]),showWarnings=F,recursive=T )
  pdf( file=paste0("PDF/",inp.mirna.dt[["N"]][ imir ],"/",OUTNAMES[ibf],"-",pepDate(),".pdf"),width=25,height=10,pointsize=15 )
  tmp.dt  <-
   cleavageEfficiency(
    CRDS   = unlist(strsplit( x=inp.mirna.dt[["C"]][ imir ],split=",",fixed=T )),
    PRELOC = inp.mirna.dt[["L"]][ imir ],
    CRDDEV = 2,
    xGR    = tmp.gr,
    PLOT   = T,
    TIT    = inp.mirna.dt[["N"]][ imir ]
   )
  dev.off()
  tmp.dt[["sample"]] <- OUTNAMES[ibf]
  cEff.dt <- rbind( cEff.dt,tmp.dt )
 }

}

save( list=c("cEff.dt"),file=paste0("cEff-",pepDate(),".dt.rda") )

} else   {

#load( "cEff-20220202.dt.rda",verbose=T )
load( "cEff-20220210.dt.rda",verbose=T )

}

cEff.mean.dt <- data.table()
for ( imir in seq_along(inp.mirna.dt[["N"]]) ) {
 cat( "  >> ",inp.mirna.dt[["N"]][ imir ],"\n",sep="" )

 cat( "   >> all samples\n",sep="" )

 pdf( file=paste0("PDF/",inp.mirna.dt[["N"]][ imir ],"/E15.5-all-samples-comparison-",pepDate(),".pdf"),width=3*25,height=2*10,pointsize=25 )
  plot_cEff(
   par.dt=cEff.dt[ name==inp.mirna.dt[["N"]][ imir ] ],
   WT=c("DicerXXE1","DicerXXE6","DicerXXE8B"),KO=c("DicerXXE10B","DicerXXE2","DicerXXE3B","DicerXXE4","DicerXXE7B"),
   MFROW=NULL,TIT=NULL,YLABS=setNames( c("wt","DicerX"),c("WT","KO") )
  )
 dev.off()
 pdf( file=paste0("PDF/",inp.mirna.dt[["N"]][ imir ],"/ESC-all-samples-comparison-",pepDate(),".pdf"),width=3*25,height=2*10,pointsize=25 )
  plot_cEff(
   par.dt=cEff.dt[ name==inp.mirna.dt[["N"]][ imir ] ],
   WT=c("RS7_Mos_r1","RS7_Mos_r2","RS7_Mos_r3"),KO=c("RS10_Mos_r1","RS10_Mos_r2","RS10_Mos_r3"),
   MFROW=NULL,TIT=NULL,YLABS=setNames( c("wt","DicerX"),c("WT","KO") )
  )
 dev.off()
 pdf( file=paste0("PDF/",inp.mirna.dt[["N"]][ imir ],"/PKR-all-samples-comparison-",pepDate(),".pdf"),width=3*25,height=2*10,pointsize=25 )
  plot_cEff(
   par.dt=cEff.dt[ name==inp.mirna.dt[["N"]][ imir ] ],
   WT=c("embryo_Prkra_WT_r1","embryo_Prkra_WT_r2"),KO=c("embryo_Prkra_Mut_r1","embryo_Prkra_Mut_r2","embryo_Prkra_Mut_r3"),
   MFROW=NULL,TIT=NULL,YLABS=setNames( c("wt","Prkra-ko"),c("WT","KO") )
  )
 dev.off()
 pdf( file=paste0("PDF/",inp.mirna.dt[["N"]][ imir ],"/TAR-all-samples-comparison-",pepDate(),".pdf"),width=3*25,height=2*10,pointsize=25 )
  plot_cEff(
   par.dt=cEff.dt[ name==inp.mirna.dt[["N"]][ imir ] ],
   WT=c("embryo_Tarbp2_WT_r1","embryo_Tarbp2_WT_r2","embryo_Tarbp2_WT_r3"),KO=c("embryo_Tarbp2_Mut_r1","embryo_Tarbp2_Mut_r2","embryo_Tarbp2_Mut_r3"),
   MFROW=NULL,TIT=NULL,YLABS=setNames( c("wt","Tarbp2-ko"),c("WT","KO") )
  )
 dev.off()
 pdf( file=paste0("PDF/",inp.mirna.dt[["N"]][ imir ],"/GNT-all-samples-comparison-",pepDate(),".pdf"),width=3*25,height=2*10,pointsize=25 )
  plot_cEff(
   par.dt=cEff.dt[ name==inp.mirna.dt[["N"]][ imir ] ],
   WT=c("x11_WT","x14_WT","x16_WT"),KO=c("x3_GNT_homozygot","x4_GNT_homozygot","x9_GNT_homozygot"),
   MFROW=NULL,TIT=NULL,YLABS=setNames( c("wt","DicerX-GNTmut"),c("WT","KO") )
  )
 dev.off()

 cat( "   >> avr samples\n",sep="" )

 pdf( file=paste0("PDF/",inp.mirna.dt[["N"]][ imir ],"/E15.5-MEAN-comparison-",       pepDate(),".pdf"),width=2*25,height=2*10,pointsize=25 )
  qqq <-
   avrPlot_cEff(
    par.dt=cEff.dt[ name==inp.mirna.dt[["N"]][ imir ] ],
    WT=c("DicerXXE1","DicerXXE6","DicerXXE8B"),KO=c("DicerXXE10B","DicerXXE2","DicerXXE3B","DicerXXE4","DicerXXE7B"),
    MFROW=c(1,2),YLABS=setNames( c("wt","DicerX"),c("WT","KO") )
   )
  qqq[["library"]] <- "E15.5"
 dev.off()
 cEff.mean.dt <- rbind( cEff.mean.dt,qqq )
 pdf( file=paste0("PDF/",inp.mirna.dt[["N"]][ imir ],"/ESC-MEAN-comparison-",       pepDate(),".pdf"),width=2*25,height=2*10,pointsize=25 )
  qqq <-
   avrPlot_cEff(
    par.dt=cEff.dt[ name==inp.mirna.dt[["N"]][ imir ] ],
    WT=c("RS7_Mos_r1","RS7_Mos_r2","RS7_Mos_r3"),KO=c("RS10_Mos_r1","RS10_Mos_r2","RS10_Mos_r3"),
    MFROW=c(1,2),YLABS=setNames( c("wt","DicerX"),c("WT","KO") )
   )
  qqq[["library"]] <- "ESC"
 dev.off()
 cEff.mean.dt <- rbind( cEff.mean.dt,qqq )
 pdf( file=paste0("PDF/",inp.mirna.dt[["N"]][ imir ],"/PKR-MEAN-comparison-",       pepDate(),".pdf"),width=2*25,height=2*10,pointsize=25 )
  qqq <-
   avrPlot_cEff(
    par.dt=cEff.dt[ name==inp.mirna.dt[["N"]][ imir ] ],
    WT=c("embryo_Prkra_WT_r1","embryo_Prkra_WT_r2"),KO=c("embryo_Prkra_Mut_r1","embryo_Prkra_Mut_r2","embryo_Prkra_Mut_r3"),
    MFROW=c(1,2),YLABS=setNames( c("wt","Prkra-ko"),c("WT","KO") )
   )
  qqq[["library"]] <- "PKR"
 dev.off()
 cEff.mean.dt <- rbind( cEff.mean.dt,qqq )
 pdf( file=paste0("PDF/",inp.mirna.dt[["N"]][ imir ],"/TAR-MEAN-comparison-",       pepDate(),".pdf"),width=2*25,height=2*10,pointsize=25 )
  qqq <-
   avrPlot_cEff(
    par.dt=cEff.dt[ name==inp.mirna.dt[["N"]][ imir ] ],
    WT=c("embryo_Tarbp2_WT_r1","embryo_Tarbp2_WT_r2","embryo_Tarbp2_WT_r3"),KO=c("embryo_Tarbp2_Mut_r1","embryo_Tarbp2_Mut_r2","embryo_Tarbp2_Mut_r3"),
    MFROW=c(1,2),YLABS=setNames( c("wt","Tarbp2-ko"),c("WT","KO") )
   )
  qqq[["library"]] <- "TAR"
 dev.off()
 cEff.mean.dt <- rbind( cEff.mean.dt,qqq )
 pdf( file=paste0("PDF/",inp.mirna.dt[["N"]][ imir ],"/GNT-MEAN-comparison-",       pepDate(),".pdf"),width=2*25,height=2*10,pointsize=25 )
  qqq <-
   avrPlot_cEff(
    par.dt=cEff.dt[ name==inp.mirna.dt[["N"]][ imir ] ],
    WT=c("x11_WT","x14_WT","x16_WT"),KO=c("x3_GNT_homozygot","x4_GNT_homozygot","x9_GNT_homozygot"),
    MFROW=c(1,2),YLABS=setNames( c("wt","DicerX-GNTmut"),c("WT","KO") )
   )
  qqq[["library"]] <- "GNT"
 dev.off()
 cEff.mean.dt <- rbind( cEff.mean.dt,qqq )

}

print( warnings() )

#> create a new open work-book
xwb <- createWorkbook()
#> add new sheet to work-book "xwb" : max 31 characters !!!
addWorksheet(   wb=xwb,sheetName="means" )
addWorksheet(   wb=xwb,sheetName="all.replicates" )
#> close work-book "xwb"
writeDataTable( wb=xwb,sheet=    "means",          x=cEff.mean.dt,   colNames=T,rowNames=F,headerStyle=createStyle(textDecoration="bold") )
writeDataTable( wb=xwb,sheet=    "all.replicates", x=cEff.dt,        colNames=T,rowNames=F,headerStyle=createStyle(textDecoration="bold") )

#> save work-book "xwb" into file
saveWorkbook(   wb=xwb,file=paste0("partial-processing-analysis-",pepDate(),".xlsx"),overwrite=T )

