RLIBLOC="/storage/brno2/home/pepap/R/x86_64-pc-linux-gnu-library/3.2"

library( package = "S4Vectors",                    lib.loc = RLIBLOC )
library( package = "GenomeInfoDb",                 lib.loc = RLIBLOC )
library( package = "GenomicRanges",                lib.loc = RLIBLOC )
library( package = "matrixStats",                  lib.loc = RLIBLOC )
library( package = "DelayedArray",                 lib.loc = RLIBLOC )
library( package = "SummarizedExperiment",         lib.loc = RLIBLOC )
library( package = "Rsamtools",                    lib.loc = RLIBLOC )
library( package = "GenomicAlignments",            lib.loc = RLIBLOC )
library( package = "rtracklayer",                  lib.loc = RLIBLOC )
library( package = "data.table",                   lib.loc = RLIBLOC )

#= BAM file absolute path
xFILE="/storage/budejovice1/home/pepap/2017-Split/smallRNAs/s_MII_GL2015.Aligned.sortedByCoord.out.bam"
#= name for outputs
LIBNAME="GL2015"
#= minimal coverage := minimal nucleotide coverage for accepting region as "covered"
MINCOV=10
#= minimal width of covered regions
MINREGLEN=10
#= minimal gap width := regions with inter-distance <= MINGAP will be merged together
MINGAP=50
#= minimal overlap between reads and covered regions
MINOVRL=10
#= minimal aligned width of any read
MINRLEN=19
#= maximal aligned width of any read
MAXRLEN=24
#= minimal NH tag value
MINNMULTI=2
#= maximal NH tag value
MAXNMULTI=100
#= maximal number of reads to be processed at once ( saves demanded memory )
MAXRNUM=1e06

#====== SCRIPT MAIN BODY =====#

#<<<<<
#= Function := converts covered regions from BAM to GRanges

covReg <- function(gr,minCov=1,minGap=0,minLen=1,stranded=T,gene_id_prefix="cr") {

 if ( stranded ) {

#= Forward-strand covered regions
  f_gr <- GRanges(coverage(gr[ as.character(strand(gr)) == "+" ]))
  f_gr <- f_gr[ f_gr@elementMetadata$score >= minCov ]
  strand(f_gr) <- "+"
  f_gr <- reduce(x = f_gr, min.gapwidth = minGap)
#= Reverse-strand covered regions
  r_gr <- GRanges(coverage(gr[ as.character(strand(gr)) == "-" ]))
  r_gr <- r_gr[ r_gr@elementMetadata$score >= minCov ]
  strand(r_gr) <- "-"
  r_gr <- reduce(x = r_gr, min.gapwidth = minGap)

#= Join regions
  a_gr <- append(x = f_gr, values = r_gr)

 } else          {

#= All-strands covered regions
  a_gr <- GRanges(coverage(gr))
  a_gr <- a_gr[ a_gr@elementMetadata$score >= minCov ]
  strand(a_gr) <- "*"
  a_gr <- reduce(x = a_gr, min.gapwidth = minGap)

 }

#= Filter out regions which have width < minLen
 if ( minLen > 1 ) {

  a_gr <- a_gr[ width(a_gr) >= minLen ]

 }

#= Order the resulted covered regions
 covReg <- a_gr[ order(as.character(seqnames(a_gr)), start(a_gr)) ]
#= Add "gene_id" name to the covered regions
 covReg@elementMetadata$gene_id <- sprintf( fmt = paste(gene_id_prefix,"%0",nchar(length(covReg)),"d",sep=""), seq(from = 1, to = length(covReg)) )

return(covReg)

}
#>>>>>

cat(" >>> Loading BAM <<<\n",sep="")
tmpBAM <-
 readGAlignments(
  file  = xFILE,
  param = ScanBamParam( tag = c("NH","NM"), what=c("qname") )
 )

cat(" >>> Loading covReg regions <<<\n",sep="")
tmpCOVREG.gr <-
 covReg(
  gr             = tmpBAM,
  minCov         = MINCOV,
  minGap         = MINGAP,
  minLen         = MINREGLEN,
  stranded       = T,
  gene_id_prefix = LIBNAME
 )
assign(  x = paste(LIBNAME,"_CovReg.gr",sep=""), value = tmpCOVREG.gr )
save( list = paste(LIBNAME,"_CovReg.gr",sep=""), file  = paste(LIBNAME,"_CovReg.gr.rda",sep="") )
#load( file = paste(LIBNAME,"_CovReg.gr.rda",sep=""), verbose = T )
#tmpCOVREG.gr <- get(paste(LIBNAME,"_CovReg.gr",sep=""))

cat(" >>> Read selection <<<\n",sep="")
tmpBAM <-
 tmpBAM[
  ( tmpBAM@elementMetadata$NH >= MINNMULTI ) &
  ( tmpBAM@elementMetadata$NH <= MAXNMULTI ) &
  (             width(tmpBAM) >= MINRLEN )   &
  (             width(tmpBAM) <= MAXRLEN )
 ]

cat(" >>> Select reads overlapping covered regions <<<\n",sep="")
xOVR  <- findOverlaps( query = tmpCOVREG.gr, subject = tmpBAM, minoverlap = MINOVRL, ignore.strand = T )

cat(" >>> BAM reads => data.table <<<\n",sep="")
tmpDT <- data.table( qname = as.character(       tmpBAM@elementMetadata$qname[   subjectHits(xOVR) ] ),
                     reg   = as.character( tmpCOVREG.gr@elementMetadata$gene_id[   queryHits(xOVR) ] ),
                     key   = "qname" )
tmpDT <- tmpDT[ , SEL := ( length(unique(reg)) > 1 ) , by = qname ]
tmpDT <- tmpDT[   SEL == T ]
tmpDT <- unique(tmpDT)
tmpDT[["SEL"]] <- NULL
assign(  x = paste(LIBNAME,"_qname2covReg.dt",sep=""), value = tmpDT )
save( list = paste(LIBNAME,"_qname2covReg.dt",sep=""), file = paste(LIBNAME,"_qname2covReg.dt.rda",sep="") )
#load( file = paste(LIBNAME,"_qname2covReg.dt.rda",sep=""), verbose = T )
#tmpDT <- get(paste(LIBNAME,"_qname2covReg.dt",sep=""))

tmpQNAME <- unique( tmpDT[ , qname ] )
if ( length(tmpQNAME) > MAXRNUM ) {

 tmpSPLIT <- seq( from = 1, to = length(tmpQNAME), length.out = floor( length(tmpQNAME) / MAXRNUM ) )[-1]
 count.dt <- data.table()
 cat(" >>> Processing ... <<<\n",sep="")

 j <- 0
 for ( i in tmpSPLIT ) {

  j <- ( j + 1 )

  cat("  >>> ",i," / ",length(tmpQNAME)," reads <<< \n",sep="")
  cat("  >>> BAM reads => data.table <<<\n",sep="")
  tmpDT.gr <-
   GRanges(
    seqnames = tmpDT[ qname %in% tmpQNAME[ seq(j,i) ] , qname ],
    ranges   = IRanges(
     start   = 1,
     end     = 1
    ),
    strand   = "*",
    reg      = tmpDT[ qname %in% tmpQNAME[ seq(j,i) ] , reg ]
   )

  cat("  >>> Joining covered regions by reads <<<\n",sep="")
  xOVR <- findOverlaps( tmpDT.gr, drop.self = T, drop.redundant = T, ignore.strand = T )
  out.dt <-
   data.table(
    cr1   = tmpDT.gr@elementMetadata$reg[   queryHits(xOVR) ],
    cr2   = tmpDT.gr@elementMetadata$reg[ subjectHits(xOVR) ],
    int1  = as.integer( substr( x = tmpDT.gr@elementMetadata$reg[   queryHits(xOVR) ], start = nchar(LIBNAME)+1, stop = 1000 ) ),
    int2  = as.integer( substr( x = tmpDT.gr@elementMetadata$reg[ subjectHits(xOVR) ], start = nchar(LIBNAME)+1, stop = 1000 ) )
   )
  out.dt <-
   rbind(
    out.dt[ int1 <= int2 , c("cr1","cr2") , with = F ],
    data.table(
     cr1  = out.dt[ int1 > int2 , cr2 ],
     cr2  = out.dt[ int1 > int2 , cr1 ]
    )
   )

 out.dt[["reg_con"]] <- paste(out.dt[,cr1],out.dt[,cr2],sep=":")
 count.dt <-
 rbind(
  count.dt,
  out.dt[ , { read_connection_count = .N; list( read_connection_count = read_connection_count ) } , by = "reg_con" ]
 )

 rm( list = c("tmpDT.gr","xOVR") )
 gc( verbose = T )
 j <- i

 }

 count.dt <- count.dt[ , { read_connection_count = sum(read_connection_count); list( read_connection_count = read_connection_count ) } , by = "reg_con" ]

} else                            {

 cat(" >>> BAM reads => data.table <<<\n",sep="")
 tmpDT.gr <-
  GRanges(
   seqnames = tmpDT[ , qname ],
   ranges   = IRanges(
    start   = 1,
    end     = 1
   ),
   strand   = "*",
   reg      = tmpDT[ , reg ]
  )
 print(length(tmpDT.gr))

 cat(" >>> Joining covered regions by reads <<<\n",sep="")
 xOVR <- findOverlaps( tmpDT.gr, drop.self = T, drop.redundant = T, ignore.strand = T )
 out.dt <-
  data.table(
   cr1   = tmpDT.gr@elementMetadata$reg[   queryHits(xOVR) ],
   cr2   = tmpDT.gr@elementMetadata$reg[ subjectHits(xOVR) ],
   int1  = as.integer( substr( x = tmpDT.gr@elementMetadata$reg[   queryHits(xOVR) ], start = nchar(LIBNAME)+1, stop = 1000 ) ),
   int2  = as.integer( substr( x = tmpDT.gr@elementMetadata$reg[ subjectHits(xOVR) ], start = nchar(LIBNAME)+1, stop = 1000 ) )
  )
 out.dt <-
  rbind(
   out.dt[ int1 <= int2  , c("cr1","cr2") , with = F ],
   data.table( cr1 = out.dt[ int1 > int2 , cr2 ], cr2 = out.dt[ int1 > int2 , cr1 ] )
  )

out.dt[["reg_con"]] <- paste(out.dt[,cr1],out.dt[,cr2],sep=":")
count.dt <- out.dt[ , { read_connection_count = .N; list( read_connection_count = read_connection_count ) } , by = "reg_con" ]

}

assign(  x = paste(LIBNAME,"_tab.readConnections.dt",sep=""), value = count.dt )
save( list = paste(LIBNAME,"_tab.readConnections.dt",sep=""), file  = paste(LIBNAME,"_tab.readConnections.dt.rda",sep="") )
#load( file = paste(LIBNAME,"_tab.readConnections.dt.rda",sep=""), verbose = T )

#====== SCRIPT MAIN BODY =====#

