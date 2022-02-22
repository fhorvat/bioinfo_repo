library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(stringr)

#>>> FUNCTIONs
gr2dt <- function(gr) {

 dt <-
     data.table( chr    = as.character(seqnames(gr)),
                 start  = start(gr),
                 end    =   end(gr),
                 strand = as.character(strand(gr)) )
 if ( length(gr@elementMetadata) != 0 ) {
  md <-
  as.data.table( gr@elementMetadata )
  dt <- cbind(dt,md)
 }

 return(dt)

}
locSeq <- function( mirSEQ,preSEQ,DELTA=15 ) {

 tmp.list <- stringr::str_locate_all( pattern=mirSEQ,string=preSEQ )
 out.str  <- rep.int( x="check_it",times=length(tmp.list) )
 tmp.half <- round( nchar(preSEQ)/2 )
 uni.res  <- unlist(lapply( X=tmp.list,FUN=nrow ))

# return( tmp.list )

 uni.beg  <- unlist(lapply( X=tmp.list[ uni.res==1 ],FUN=function(M){ return(M[1,1]) } ))
 uni.end  <- unlist(lapply( X=tmp.list[ uni.res==1 ],FUN=function(M){ return(M[1,2]) } ))

 out.str[ uni.res==1 ][ (uni.beg+DELTA) <= tmp.half[ uni.res==1 ] ] <- "5p"
 out.str[ uni.res==1 ][ (uni.end-DELTA) >= tmp.half[ uni.res==1 ] ] <- "3p"

 return( list( BEG=uni.beg,END=uni.end,LEN=tmp.half,CNT=uni.res,STR=out.str ) )

}
#<<<

#>>> DOWLOAD INPUT FILES
system( command="wget https://www.mirbase.org/ftp/CURRENT/genomes/mmu.gff3" )
system( command="wget https://www.mirbase.org/ftp/CURRENT/mature.fa.gz"     )
system( command="wget https://www.mirbase.org/ftp/CURRENT/aliases.txt.gz"   )
#<<<

mmu.gr <- import.gff3( con = "mmu.gff3" )

if ( sum( !( levels(mmu.gr@elementMetadata$type) %in% c("miRNA_primary_transcript","miRNA") ) ) != 0 ) {
 warning("Unexpected \"type\" factors !")
}

MIRs <- mmu.gr[ as.character(mmu.gr@elementMetadata$type) == "miRNA" ]
PREs <- mmu.gr[ as.character(mmu.gr@elementMetadata$type) == "miRNA_primary_transcript" ]
xOVR <- findOverlaps( query = MIRs, subject = PREs, type = "within", ignore.strand = F )

if ( sum( duplicated(queryHits(xOVR)) ) != 0 ) {
 warning("Some mature miRNAs overlap with multiple pre-mature miRNAs !")
}

if ( length(queryHits(xOVR)) != length(MIRs) ) {
 warning("Some mature miRNAs do not align within any pre-mature miRNAs !")
}

mature_miRNAs_mmu.gr                 <- granges(MIRs)
mcols(mature_miRNAs_mmu.gr)["ID"]    <- MIRs@elementMetadata[,"ID"]
mcols(mature_miRNAs_mmu.gr)["preID"] <- MIRs@elementMetadata[,"Derives_from"]
mcols(mature_miRNAs_mmu.gr)["Name"]  <- MIRs@elementMetadata[,"Name"]

mirAnnot.dt           <- gr2dt(mature_miRNAs_mmu.gr)
mirAnnot.dt[["gLoc"]] <- mirAnnot.dt[,paste(chr,":",start,"-",end,sep="")]

mirAnnot.dt[["mirAnnotNum"]]                                       <- gsub("^mmu[-]miR[-]","",mirAnnot.dt[,Name])
mirAnnot.dt[ grepl("^mmu[-]let[-]",mirAnnotNum) ][["mirAnnotNum"]] <- as.character("0")
mirAnnot.dt[["mirAnnotNum"]]                                       <- gsub("[a-z]","", mirAnnot.dt[,mirAnnotNum])
mirAnnot.dt[["mirAnnotNum"]]                                       <- gsub("[A-Z]","", mirAnnot.dt[,mirAnnotNum])
mirAnnot.dt[["mirAnnotNum"]]                                       <- gsub("[-].*$","",mirAnnot.dt[,mirAnnotNum])
mirAnnot.dt[["mirAnnotNum"]]                                       <- as.numeric(mirAnnot.dt[["mirAnnotNum"]])
mirAnnot.dt[["tmp.ID"]] <- gsub( "_[0-9]*$","",mirAnnot.dt[["ID"]] )

aliases.dt  <- fread( "aliases.txt.gz",header=F )
mirAnnot.dt <- merge(mirAnnot.dt,aliases.dt,by.x="tmp.ID",by.y="V1",all.x=T,sort=F)
mirAnnot.dt[["miRNA"]]  <- gsub( "[-][3,5]p$","",mirAnnot.dt[["Name"]] )
mirAnnot.dt[["type"]]   <- as.character("undetermined")
mirAnnot.dt[ mirAnnot.dt[ , { tmp.name <- sub( "[-][3,5]p$","",Name ); list( STAR=grepl( pattern=paste(tmp.name,";",   sep=""),x=V2 ) ) } , by="ID" ][ , STAR ] ][["type"]] <- "miRNA"
mirAnnot.dt[ mirAnnot.dt[ , { tmp.name <- sub( "[-][3,5]p$","",Name ); list( STAR=grepl( pattern=paste(tmp.name,"[*];",sep=""),x=V2 ) ) } , by="ID" ][ , STAR ] ][["type"]] <- "miRNA*"
mirAnnot.dt[ mirAnnot.dt[ , {                                          list( STAR=grepl( pattern=paste(    Name,";",   sep=""),x=V2 ) ) } , by="ID" ][ , STAR ] & type=="undetermined" ][["type"]] <-
 "miRNA"
if ( nrow(mirAnnot.dt[ mirAnnot.dt[ , {                                list( STAR=grepl( pattern=paste(    Name,"[*];",sep=""),x=V2 ) ) } , by="ID" ][ , STAR ] & type=="undetermined" ])!=0 ) {
          mirAnnot.dt[ mirAnnot.dt[ , {                                list( STAR=grepl( pattern=paste(    Name,"[*];",sep=""),x=V2 ) ) } , by="ID" ][ , STAR ] & type=="undetermined" ][["type"]] <-
 "miRNA*"
}

miSEQ       <- readRNAStringSet("mature.fa.gz")
miSEQ       <- data.table( tmp.ID=unlist( lapply( X=names(miSEQ),FUN=function(L){ return( unlist(strsplit(x=L,split=" ",fixed=T))[2] ) } ) ),SEQ=as.character(miSEQ) )
mirAnnot.dt <- merge( mirAnnot.dt,miSEQ,by="tmp.ID",all.x=T,sort=F )

preID.crd                   <- mirAnnot.dt[ , { list( preID.gLoc=paste( unique(chr),":",min(start),"-",max(end),":",unique(strand),sep="" ) ) } , by="preID" ]
preID.seq                   <- getSeq( Mmusculus,GRanges(preID.crd[,preID.gLoc]) )
preID.crd[["pre.SEQ"]]      <- gsub( "T","U",as.character(preID.seq) )
mirAnnot.dt                 <- merge( mirAnnot.dt,preID.crd,by="preID",all.x=T,sort=F )
mirAnnot.dt[["preID.gLoc"]] <- gsub( "[:][+,-]$","",mirAnnot.dt[["preID.gLoc"]] )

mirAnnot.dt[["miRNA.STR"]]                         <- as.character("unknown")
mirAnnot.dt[ grepl("[-]3p$",Name) ][["miRNA.STR"]] <- "3p"
mirAnnot.dt[ grepl("[-]5p$",Name) ][["miRNA.STR"]] <- "5p"

mirAnnot.dt[["V2"]]     <- NULL
mirAnnot.dt[["tmp.ID"]] <- NULL

mir.str.list <- locSeq( mirSEQ=mmu.all.dt[["SEQ"]],preSEQ=mmu.all.dt[["preSEQ"]],DELTA=15 )
mmu.all.dt[["det.STR"]] <- mir.str.list$STR

mmu.all.dt[ ID=="MIMAT0020633" ][["det.STR"]] <- "3p"
mmu.all.dt[ ID=="MIMAT0017329" ][["det.STR"]] <- "5p"
mmu.all.dt[ ID=="MIMAT0017325" ][["det.STR"]] <- "5p"

print( table(mmu.all.dt[["det.STR"]]) )

mirAnnot.dt <- merge( mirAnnot.dt,mmu.all.dt[,c("ID","det.STR"),with=F],by="ID",all=T,sort=F )

PEPFAC=2
MINEXP=0.1
CONSIDERED_EXPERIMENTS=c( "pep.type.EMB","pep.type.GNT" )

#= EMB15.5-ES
load( paste0("EMB15.5/CSVs",SHROF,".rda"),verbose=T )
emb15.5.RPM.dt <- data.table( ID=ddd[["ID"]],WT.means=log10(ddd[["baseMean"]]+1) )
rm( list=c("rrr","ddd") )

mirAnnot.dt[["pep.type.EMB"]] <- as.character("notClear")
undef.dt                      <- premirAnnot.dt
xtype.dt                      <- merge( undef.dt,emb15.5.RPM.dt,by.x="ID.5p",by.y="ID",all.x=T,sort=F                         )
xtype.dt                      <- merge( xtype.dt,emb15.5.RPM.dt,by.x="ID.3p",by.y="ID",all.x=T,sort=F,suffixes=c(".5p",".3p") )
mirAnnot.dt[ ID %in% xtype.dt[ ( WT.means.5p *   PEPFAC ) < WT.means.3p  , ID.5p ] ][["pep.type.EMB"]] <- "miRNA*"
mirAnnot.dt[ ID %in% xtype.dt[   WT.means.5p > ( PEPFAC *   WT.means.3p ), ID.3p ] ][["pep.type.EMB"]] <- "miRNA*"
mirAnnot.dt[ ID %in% xtype.dt[ ( WT.means.5p *   PEPFAC ) < WT.means.3p  , ID.3p ] ][["pep.type.EMB"]] <- "miRNA"
mirAnnot.dt[ ID %in% xtype.dt[   WT.means.5p > ( PEPFAC *   WT.means.3p ), ID.5p ] ][["pep.type.EMB"]] <- "miRNA"
mirAnnot.dt[ ID %in% xtype.dt[ ( WT.means.5p<MINEXP ) & ( WT.means.3p<MINEXP ) , ID.5p ] ][["pep.type.EMB"]] <- "lowExp"
mirAnnot.dt[ ID %in% xtype.dt[ ( WT.means.5p<MINEXP ) & ( WT.means.3p<MINEXP ) , ID.3p ] ][["pep.type.EMB"]] <- "lowExp"

#= EMB15.5-DcrGNT-VB
load( paste0("GNT/CSVs",SHROF,".rda"),verbose=T )
gnt.RPM.dt <- data.table( ID=ddd[["ID"]],WT.means=log10(ddd[["baseMean"]]+1) )
rm( list=c("rrr","ddd") )

mirAnnot.dt[["pep.type.GNT"]] <- as.character("notClear")
undef.dt                      <- premirAnnot.dt
xtype.dt                      <- merge( undef.dt,gnt.RPM.dt,by.x="ID.5p",by.y="ID",all.x=T,sort=F                         )
xtype.dt                      <- merge( xtype.dt,gnt.RPM.dt,by.x="ID.3p",by.y="ID",all.x=T,sort=F,suffixes=c(".5p",".3p") )
mirAnnot.dt[ ID %in% xtype.dt[ ( WT.means.5p *   PEPFAC ) < WT.means.3p  , ID.5p ] ][["pep.type.GNT"]] <- "miRNA*"
mirAnnot.dt[ ID %in% xtype.dt[   WT.means.5p > ( PEPFAC *   WT.means.3p ), ID.3p ] ][["pep.type.GNT"]] <- "miRNA*"
mirAnnot.dt[ ID %in% xtype.dt[ ( WT.means.5p *   PEPFAC ) < WT.means.3p  , ID.3p ] ][["pep.type.GNT"]] <- "miRNA"
mirAnnot.dt[ ID %in% xtype.dt[   WT.means.5p > ( PEPFAC *   WT.means.3p ), ID.5p ] ][["pep.type.GNT"]] <- "miRNA"
mirAnnot.dt[ ID %in% xtype.dt[ ( WT.means.5p<MINEXP ) & ( WT.means.3p<MINEXP ) , ID.5p ] ][["pep.type.GNT"]] <- "lowExp"
mirAnnot.dt[ ID %in% xtype.dt[ ( WT.means.5p<MINEXP ) & ( WT.means.3p<MINEXP ) , ID.3p ] ][["pep.type.GNT"]] <- "lowExp"

mirAnnot.pepType.dt <- mirAnnot.dt

#= (1) miRNA/miRNA* from miRbase

undef.dt                  <- premirAnnot.dt[ ( type.5p==type.3p ) & ( ( type.3p=="miRNA" ) | ( type.3p=="miRNA*" ) ) ]

MIRBASE <-
 c(
  premirAnnot.dt[ ( ( type.5p=="miRNA" ) & ( type.3p=="miRNA*" ) ) | ( ( type.5p=="miRNA*" ) & ( type.3p=="miRNA" ) ) , ID.5p ],
  premirAnnot.dt[ ( ( type.5p=="miRNA" ) & ( type.3p=="miRNA*" ) ) | ( ( type.5p=="miRNA*" ) & ( type.3p=="miRNA" ) ) , ID.3p ]
 )
XSTRAND <- c( undef.dt[,ID.5p],undef.dt[,ID.3p] )
ONLYONE <- mirAnnot.pepType.dt[ preID %in% names(table(mirAnnot.pepType.dt[["preID"]]))[ table(mirAnnot.pepType.dt[["preID"]])==1 ] ][["ID"]]

if ( nrow(mirAnnot.pepType.dt) != length(MIRBASE) + length(XSTRAND) + length(ONLYONE) ) { stop("\n !!! Incomplete categories !!!\n\n") }

mirAnnot.pepType.dt[["pep.type"]]                    <- mirAnnot.pepType.dt[["type"]]

print( mirAnnot.pepType.dt[ , sort(table(pep.type),decreasing=T) ] )

mirAnnot.pepType.dt[ ID %in% ONLYONE ][["pep.type"]] <- as.character("single_miRNA")
mirAnnot.pepType.dt[ ID %in% XSTRAND ][["pep.type"]] <-
 apply(
  X      = as.matrix( x=mirAnnot.pepType.dt[ ID %in% XSTRAND , CONSIDERED_EXPERIMENTS , with=F ] ),
  MARGIN = 1,
  FUN    = function(LINE) {
   return( paste(sort(unique(LINE)),collapse="|") )
  }
 )

PEP.ANNOT <-
 setNames(
  object = c(
         "miRNA",        "miRNA*",       "notClear",        "notClear", "cellSpecific",
         "cellSpecific",          "cellSpecific",                 "cellSpecific",
               "notClear",               "notClear",        "notClear"
            ),
  nm     = c(
  "lowExp|miRNA", "lowExp|miRNA*", "miRNA|notClear", "miRNA*|notClear", "miRNA|miRNA*",
  "lowExp|miRNA|miRNA*", "miRNA|miRNA*|notClear", "lowExp|miRNA|miRNA*|notClear",
  "lowExp|miRNA|notClear", "lowExp|miRNA*|notClear", "lowExp|notClear"
            )
 )

for ( icase in seq_along(PEP.ANNOT) ) {
 if ( nrow(mirAnnot.pepType.dt[ pep.type==names(PEP.ANNOT)[icase] ])!=0 ) {
  mirAnnot.pepType.dt[ pep.type==names(PEP.ANNOT)[icase] ][["pep.type"]] <- as.character( as.vector(PEP.ANNOT)[icase] )
 }
}

print( mirAnnot.pepType.dt[ , sort(table(pep.type),decreasing=T) ] )

MIRTRONS                                                               <- fread( "mirtrons-mirBase.csv",header=T )
mirAnnot.dt                                                            <- mirAnnot.pepType.dt
mirAnnot.dt[["mirtron"]]                                               <- as.logical("FALSE")
mirAnnot.dt[ preID %in% unique(MIRTRONS[["Accession"]]) ][["mirtron"]] <- as.logical("TRUE")

save( mirAnnot.dt,file="mirAnnot.dt.rda" )

