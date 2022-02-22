library(data.table)
library(GenomicAlignments)

load( "mirAnnot.dt.rda",verbose=T )

mirAnnot.dt[["CleavagePoint"]]                  <- as.numeric("")
mirAnnot.dt[  det.STR=="5p" & strand=="+" ][["CleavagePoint"]] <-
 mirAnnot.dt[ det.STR=="5p" & strand=="+" ,   end ]
mirAnnot.dt[  det.STR=="3p" & strand=="+" ][["CleavagePoint"]] <-
 mirAnnot.dt[ det.STR=="3p" & strand=="+" , start ]
mirAnnot.dt[  det.STR=="3p" & strand=="-" ][["CleavagePoint"]] <-
 mirAnnot.dt[ det.STR=="3p" & strand=="-" ,   end ]
mirAnnot.dt[  det.STR=="5p" & strand=="-" ][["CleavagePoint"]] <-
 mirAnnot.dt[ det.STR=="5p" & strand=="-" , start ]

save( mirAnnot.dt,file="mirAnnot-CPs.dt.rda" )

RLRAN=c(19,25)
ONAME=paste0(RLRAN[1],"to",RLRAN[2],"nt")

DELTA=15

CPs.gr        <- GRanges( mirAnnot.dt[ , paste(chr,":",CleavagePoint-DELTA,"-",CleavagePoint+DELTA,":",strand,sep="") ], ID=mirAnnot.dt[ , ID ] )
CPs.5p.fwd.gr <- CPs.gr[ as.character(strand(CPs.gr))=="+" & ( mcols(CPs.gr)[["ID"]] %in% mirAnnot.dt[ det.STR=="5p", ID ] ) ]
CPs.5p.rev.gr <- CPs.gr[ as.character(strand(CPs.gr))=="-" & ( mcols(CPs.gr)[["ID"]] %in% mirAnnot.dt[ det.STR=="5p", ID ] ) ]
CPs.3p.fwd.gr <- CPs.gr[ as.character(strand(CPs.gr))=="+" & ( mcols(CPs.gr)[["ID"]] %in% mirAnnot.dt[ det.STR=="3p", ID ] ) ]
CPs.3p.rev.gr <- CPs.gr[ as.character(strand(CPs.gr))=="-" & ( mcols(CPs.gr)[["ID"]] %in% mirAnnot.dt[ det.STR=="3p", ID ] ) ]

EMB.ES.BAM <- c( "DicerX_KO_10B_r2.bam","DicerX_WT_1_r3.bam","DicerX_KO_2_r5.bam","DicerX_KO_3B_r3.bam","DicerX_KO_4_r1.bam","DicerX_WT_6_r1.bam","DicerX_KO_7B_r4.bam","DicerX_WT_8B_r2.bam" )
EMB.ES.CON <- c( "DicerX_KO_10B_r2    ","DicerX_WT_1_r3    ","DicerX_KO_2_r5    ","DicerX_KO_3B_r3    ","DicerX_KO_4_r1    ","DicerX_WT_6_r1    ","DicerX_KO_7B_r4    ","DicerX_WT_8B_r2    " )
ESC.ES.BAM <- c( "RS10_Mos_r1.bam", "RS10_Mos_r2.bam", "RS10_Mos_r3.bam", "RS7_Mos_r1.bam", "RS7_Mos_r2.bam", "RS7_Mos_r3.bam" )
ESC.ES.CON <- c( "RS10_Mos_r1",     "RS10_Mos_r2",     "RS10_Mos_r3",     "RS7_Mos_r1",     "RS7_Mos_r2",     "RS7_Mos_r3"     )
EMB.PU.BAM <- c(
 "emb_Prkra_Mut_r1.bam", "emb_Prkra_Mut_r2.bam", "emb_Prkra_Mut_r3.bam",
 "emb_Prkra_WT_r1.bam",  "emb_Prkra_WT_r2.bam",
 "emb_Tarbp2_Mut_r1.bam","emb_Tarbp2_Mut_r2.bam","emb_Tarbp2_Mut_r3.bam",
 "emb_Tarbp2_WT_r1.bam", "emb_Tarbp2_WT_r2.bam", "emb_Tarbp2_WT_r3.bam"
 )
EMB.PU.CON <- c(
 "emb_Prkra_Mut_r1", "emb_Prkra_Mut_r2", "emb_Prkra_Mut_r3",
 "emb_Prkra_WT_r1",  "emb_Prkra_WT_r2",
 "emb_Tarbp2_Mut_r1","emb_Tarbp2_Mut_r2","emb_Tarbp2_Mut_r3",
 "emb_Tarbp2_WT_r1", "emb_Tarbp2_WT_r2", "emb_Tarbp2_WT_r3"
 )
GNT.VB.BAM <- c( "x11_WT.bam","x14_WT.bam","x16_WT.bam","x3_GNT_homozygot.bam","x4_GNT_homozygot.bam","x9_GNT_homozygot.bam" )
GNT.VB.CON <- c( "x11_WT",    "x14_WT",    "x16_WT",    "x3_GNT_homozygot",    "x4_GNT_homozygot",    "x9_GNT_homozygot"     )

BAMFILES <- c( EMB.ES.BAM,ESC.ES.BAM,EMB.PU.BAM,GNT.VB.BAM )
BAMNAMES <- c( EMB.ES.CON,ESC.ES.CON,EMB.PU.CON,GNT.VB.CON )

cat("\n",sep="")
libSize.vec <- c()
out.objs    <- c()
for ( i in seq_along(BAMNAMES) ) {

 cat( " >> ", BAMNAMES[i]," : <<\n",sep="" )

print("check-01")
 tmp.ga     <- readGAlignments( file=BAMFILES[i],param=ScanBamParam( tag=c("HI","NH","nM") ) )
 tmp.ga     <- tmp.ga[ njunc(tmp.ga)==0 ]
print("check-02")
 tmp.gr               <- as( tmp.ga,"GRanges" )
 tmp.gr               <- tmp.gr[ ( width(tmp.gr) >= RLRAN[1] ) & ( width(tmp.gr) <= RLRAN[2] ) ]

 libSize.vec          <- append( libSize.vec,sum( mcols(tmp.gr)[["HI"]]==1 ) )

 tmp.5p.fwd.gr        <- tmp.gr[ as.character(strand(tmp.gr))=="+" ]
 start(tmp.5p.fwd.gr) <-   end(tmp.5p.fwd.gr)
 tmp.5p.rev.gr        <- tmp.gr[ as.character(strand(tmp.gr))=="-" ]
   end(tmp.5p.rev.gr) <- start(tmp.5p.rev.gr)
 tmp.3p.fwd.gr        <- tmp.gr[ as.character(strand(tmp.gr))=="+" ]
   end(tmp.3p.fwd.gr) <- start(tmp.3p.fwd.gr)
 tmp.3p.rev.gr        <- tmp.gr[ as.character(strand(tmp.gr))=="-" ]
 start(tmp.3p.rev.gr) <-   end(tmp.3p.rev.gr)
print("check-03")
 tmp.5p.fwd.cv <- coverage( tmp.5p.fwd.gr )
 tmp.5p.rev.cv <- coverage( tmp.5p.rev.gr )
 tmp.3p.fwd.cv <- coverage( tmp.3p.fwd.gr )
 tmp.3p.rev.cv <- coverage( tmp.3p.rev.gr )
print("check-04")
 tmp.5p.fwd.cv <- tmp.5p.fwd.cv[ CPs.5p.fwd.gr ]
 tmp.5p.rev.cv <- tmp.5p.rev.cv[ CPs.5p.rev.gr ]
 tmp.3p.fwd.cv <- tmp.3p.fwd.cv[ CPs.3p.fwd.gr ]
 tmp.3p.rev.cv <- tmp.3p.rev.cv[ CPs.3p.rev.gr ]
print("check-05")
 tmp.dt     <- rbind( as.data.table( tmp.5p.fwd.cv ),as.data.table( tmp.5p.rev.cv ),as.data.table( tmp.3p.fwd.cv ),as.data.table( tmp.3p.rev.cv ) )

print("check-06")
 tmp.mat <- matrix(
  data=tmp.dt[,value],ncol=(2*DELTA)+1,byrow=T,
  dimnames=list( c(mcols(CPs.5p.fwd.gr)[["ID"]],mcols(CPs.5p.rev.gr)[["ID"]],mcols(CPs.3p.fwd.gr)[["ID"]],mcols(CPs.3p.rev.gr)[["ID"]]),c() )
 )
 out.mat <- tmp.mat[ mirAnnot.dt[ strand=="+" , ID ] , ]
 tmp.mat <- tmp.mat[ mirAnnot.dt[ strand=="-" , ID ] , rev(seq(ncol(tmp.mat))) ]
 out.mat <- rbind( out.mat,tmp.mat )
 colnames(out.mat) <- seq( from=(-1)*DELTA,to=(+1)*DELTA )

print("check-07")
 assign( x=paste( BAMNAMES[i],".CPs.mat",sep="" ),value=out.mat )
 out.objs    <- append( out.objs,paste( BAMNAMES[i],".CPs.mat",sep="" ) )

cat("\n",sep="")

}
cat("\n",sep="")

tmp.dt <- data.table( LIBS=BAMNAMES,MAPPED=libSize.vec )
assign(  x=paste("DicerX.",ONAME,".libSize.dt",sep=""),value=tmp.dt )
save( list=paste("DicerX.",ONAME,".libSize.dt",sep=""),file=paste("DicerX.",ONAME,".libSize.dt.rda",sep="") )

save( list=out.objs,file=paste("DicerX-CPs.",ONAME,"-",pepDate(),".rda",sep="") )

