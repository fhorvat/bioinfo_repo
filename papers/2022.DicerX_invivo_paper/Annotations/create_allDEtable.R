library(data.table)

load( "mirAnnot.dt.rda",verbose=T )

SHROF    <- c(".ORIGINAL-LFC",".SHRUNKEN-LFC")[2]
COLNAMES <- c("ID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")

#= EMB15.5-ES
DEDTFILE <-
 system(
  command=paste0("/bin/ls /storage/brno1-cerit/home/pepap/Taborska_smallRNAlibs_DicerXXE_191217/BAM/00.stats/01.miRNAcounts/01.final/EMB15.5-ES-smallRNAseq",SHROF,".de.dt.rda"),
  intern=T
 )
emb15.5.DE.dt <- load( DEDTFILE,verbose=T )
emb15.5.DE.dt <- get(emb15.5.DE.dt)
emb15.5.DE.dt <- emb15.5.DE.dt[ ,COLNAMES,with=F ]
colnames(emb15.5.DE.dt)[-1] <- paste0(COLNAMES[-1],".EMB")

#= ESC-ES
DEDTFILE <-
 system(
  command=paste0("/bin/ls /storage/brno1-cerit/home/pepap/Eliska.Taborska/BAM/00.stats/01.miRNAcounts/01.final/ESC-ES-smallRNAseq",SHROF,".de.dt.rda"),
  intern=T
 )
esc.DE.dt <- load( DEDTFILE,verbose=T )
esc.DE.dt <- get(esc.DE.dt)
esc.DE.dt <- esc.DE.dt[ ,COLNAMES,with=F ]
colnames(esc.DE.dt)[-1] <- paste0(COLNAMES[-1],".ESC")

#= EMB15.5-Pu
DEDTFILE <-
 system(
  command=paste0("/bin/ls /storage/brno1-cerit/home/pepap/Pullagura_2018_Genetics_PRJNA423238/BAM/READs/00.stats/01.miRNAcounts/01.final/EMB15.5-Pullagura-smallRNAseq.PKR",SHROF,".de.dt.rda"),
  intern=T
 )
pkr.DE.dt <- load( DEDTFILE,verbose=T )
pkr.DE.dt <- get(pkr.DE.dt)
pkr.DE.dt <- pkr.DE.dt[ ,COLNAMES,with=F ]
colnames(pkr.DE.dt)[-1] <- paste0(COLNAMES[-1],".PKR")
DEDTFILE <-
 system(
  command=paste0("/bin/ls /storage/brno1-cerit/home/pepap/Pullagura_2018_Genetics_PRJNA423238/BAM/READs/00.stats/01.miRNAcounts/01.final/EMB15.5-Pullagura-smallRNAseq.TAR",SHROF,".de.dt.rda"),
  intern=T
 )
tar.DE.dt <- load( DEDTFILE,verbose=T )
tar.DE.dt <- get(tar.DE.dt)
tar.DE.dt <- tar.DE.dt[ ,COLNAMES,with=F ]
colnames(tar.DE.dt)[-1] <- paste0(COLNAMES[-1],".TAR")

#= EMB15.5-DcrGNT-VB
DEDTFILE <-
 system(
  command=paste0("/bin/ls /storage/brno1-cerit/home/pepap/Valeria.Buccheri/01.Dicer-helicase-domain-mutations/BAM.MERGED/00.stats/01.miRNAcounts/no.x19_KO/01.final/EMB15.5-DcrGNT-VB",SHROF,".de.dt.rda"),
  intern=T
 )
gnt.DE.dt <- load( DEDTFILE,verbose=T )
gnt.DE.dt <- get(gnt.DE.dt)
gnt.DE.dt <- gnt.DE.dt[ ,COLNAMES,with=F ]
colnames(gnt.DE.dt)[-1] <- paste0(COLNAMES[-1],".GNT")

all.DE.dt <- merge( emb15.5.DE.dt,esc.DE.dt,by="ID",sort=F,all=T )
all.DE.dt <- merge( all.DE.dt,    tar.DE.dt,by="ID",sort=F,all=T )
all.DE.dt <- merge( all.DE.dt,    pkr.DE.dt,by="ID",sort=F,all=T )
all.DE.dt <- merge( all.DE.dt,    gnt.DE.dt,by="ID",sort=F,all=T )

save( all.DE.dt,file=paste0("all.DE",SHROF,".dt.rda") )

