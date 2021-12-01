### INFO: 
### DATE: Tue May 29 14:56:37 2018
### AUTHOR: Pepa
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# mouse
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/mouse/connected_regions")
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# # cow
# setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/cow/connected_regions")
# genome_path <- "/common/DB/genome_reference/cow/bosTau8.UMD3.1.GCA_000003055.4"
# 
# # pig
# setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/pig/connected_regions")
# genome_path <- "/common/DB/genome_reference/pig/susScr11.Sscrofa11.1.GCA_000003025.6"

######################################################## LIBRARIES
library("S4Vectors")
library("GenomeInfoDb")
library("GenomicRanges")
library("matrixStats")
library("DelayedArray")
library("SummarizedExperiment")
library("Rsamtools")
library("GenomicAlignments")
library("rtracklayer")
library("data.table")
library(magrittr)
library(stringr)

######################################################## SOURCE FILES
# lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
# source(file.path(lib_path, "wideScreen.R"))
# source(file.path(lib_path, "headt.R"))
# source(file.path(lib_path, "asdf.R"))
# wideScreen()

######################################################## FUNCTIONS
# converts covered regions from BAM to GRanges
covReg <- function(gr, minCov = 1, minGap = 0, minLen = 1, stranded = T, gene_id_prefix = "cr"){

  if(stranded){
    
    # Forward - strand covered regions
    f_gr <- GRanges(coverage(gr[as.character(strand(gr)) == "+"]))
    f_gr <- f_gr[f_gr@elementMetadata$score >= minCov]
    strand(f_gr) <- "+"
    f_gr <- reduce(x = f_gr, min.gapwidth = minGap)
    
    # Reverse-strand covered regions
    r_gr <- GRanges(coverage(gr[ as.character(strand(gr)) == "-" ]))
    r_gr <- r_gr[r_gr@elementMetadata$score >= minCov]
    strand(r_gr) <- "-"
    r_gr <- reduce(x = r_gr, min.gapwidth = minGap)
    
    # Join regions
    a_gr <- append(x = f_gr, values = r_gr)
    
  }else{
    
    # All-strands covered regions
    a_gr <- GRanges(coverage(gr))
    a_gr <- a_gr[a_gr@elementMetadata$score >= minCov]
    strand(a_gr) <- "*"
    a_gr <- reduce(x = a_gr, min.gapwidth = minGap)
    
  }
  
  # Filter out regions which have width < minLen
  if(minLen > 1){
    a_gr <- a_gr[width(a_gr) >= minLen]
  }
  
  # Order the resulted covered regions
  covReg <- a_gr[order(as.character(seqnames(a_gr)), start(a_gr)) ]
  
  # Add "gene_id" name to the covered regions
  covReg@elementMetadata$gene_id <- sprintf(fmt = paste(gene_id_prefix, "%0", nchar(length(covReg)), "d", sep = ""), seq(from = 1, to = length(covReg)))
  
  return(covReg)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# BAM file absolute path
# xFILE <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/data_sets/mouse/Tam_2008_Nature_GSE10364/Filtered/s_oocyte_19to24.SE.perfect.bam"
xFILE <- "%BAM_PATH"

# name for outputs
LIBNAME <- basename(xFILE) %>% str_remove(., ".SE.perfect.bam")
experiment_name <- xFILE %>% str_remove(., "/Filtered/.*") %>% basename(.) %>% str_remove(., "(?<=[0-9])_.*")

# minimal coverage = minimal nucleotide coverage for accepting region as "covered"
MINCOV <- 10

# minimal width of covered regions
MINREGLEN <- 10

# minimal gap width = regions with inter-distance <= MINGAP will be merged together
MINGAP <- 50

# minimal overlap between reads and covered regions
MINOVRL <- 10

# minimal aligned width of any read
MINRLEN <- 19

# maximal aligned width of any read
MAXRLEN <- 24

# minimal NH tag value
MINNMULTI <- 2

# maximal NH tag value
MAXNMULTI <- 100

# maximal number of reads to be processed at once (saves demanded memory)
MAXRNUM <- 1e06

# minimal number of reads shared between two regions
MINRNUM <- 10

######################################################## READ DATA
# load bam
tmpBAM <- readGAlignments(file  = xFILE, param = ScanBamParam(tag = c("NH", "NM"), what = c("qname")))

######################################################## MAIN CODE
# get GRanges of covered regions
tmpCOVREG.gr <- covReg(gr = tmpBAM, minCov = MINCOV, minGap = MINGAP, minLen = MINREGLEN, stranded = T, gene_id_prefix = LIBNAME)
tmpCOVREG.dt <- as.data.table(tmpCOVREG.gr)

# select reads
tmpBAM <- tmpBAM[(tmpBAM@elementMetadata$NH >= MINNMULTI) & 
                   (tmpBAM@elementMetadata$NH <= MAXNMULTI) & 
                   (width(tmpBAM) >= MINRLEN) & 
                   (width(tmpBAM) <= MAXRLEN)]

# get reads overlaping covered regions
xOVR  <- findOverlaps(query = tmpCOVREG.gr, subject = tmpBAM, minoverlap = MINOVRL, ignore.strand = T)

# get data.table of reads 
tmpDT <- data.table(qname = as.character(tmpBAM@elementMetadata$qname[subjectHits(xOVR)]),
                    reg = as.character(tmpCOVREG.gr@elementMetadata$gene_id[queryHits(xOVR)]),
                    key = "qname")
tmpDT <- tmpDT[ , SEL := (length(unique(reg)) > 1), by = qname ]
tmpDT <- tmpDT[SEL == T]
tmpDT <- unique(tmpDT)
tmpDT[["SEL"]] <- NULL


### join regions by reads - in segments to lower memory load
tmpQNAME <- unique(tmpDT[ , qname])

if(length(tmpQNAME) > MAXRNUM){
  
  tmpSPLIT <- seq(from = 1, to = length(tmpQNAME), length.out = floor(length(tmpQNAME) / MAXRNUM))[-1]
  count.dt <- data.table()

  j <- 0
  
  for (i in tmpSPLIT) {
    
    j <- (j + 1)
    
    # prepare GRanges
    tmpDT.gr <- GRanges(seqnames = tmpDT[qname %in% tmpQNAME[seq(j, i)], qname],
                        ranges = IRanges(start = 1, end = 1),
                        strand = "*",
                        reg = tmpDT[qname %in% tmpQNAME[seq(j, i)], reg])
    
    # join covered regions
    xOVR <- findOverlaps(tmpDT.gr, drop.self = T, drop.redundant = T, ignore.strand = T)
    out.dt <- data.table(cr1 = tmpDT.gr@elementMetadata$reg[queryHits(xOVR)],
                         cr2 = tmpDT.gr@elementMetadata$reg[subjectHits(xOVR)],
                         int1 = as.integer(substr(x = tmpDT.gr@elementMetadata$reg[queryHits(xOVR)], start = nchar(LIBNAME) + 1, stop = 1000)),
                         int2 = as.integer(substr(x = tmpDT.gr@elementMetadata$reg[subjectHits(xOVR)], start = nchar(LIBNAME) + 1, stop = 1000)))
    out.dt <- rbind(out.dt[int1 <= int2, c("cr1", "cr2") , with = F],
                    data.table(cr1 = out.dt[int1 > int2, cr2],
                               cr2 = out.dt[int1 > int2, cr1]))
    
    out.dt[["reg_con"]] <- paste(out.dt[, cr1], out.dt[, cr2], sep = ":")
    
    count.dt <- rbind(count.dt, out.dt[ , {read_connection_count = .N; list(read_connection_count = read_connection_count)}, by = "reg_con"])
    
    rm(list = c("tmpDT.gr","xOVR") )
    gc(verbose = T)
    
    j <- i
    
  }
  
  count.dt <- count.dt[ , {read_connection_count = sum(read_connection_count); list(read_connection_count = read_connection_count)} , by = "reg_con"]
  
  
} else {
  
  tmpDT.gr <- GRanges(seqnames = tmpDT[ , qname],
                      ranges = IRanges(start = 1, end = 1),
                      strand = "*",
                      reg = tmpDT[ , reg])
  print(length(tmpDT.gr))
  
  # join covered regions by reads
  xOVR <- findOverlaps(tmpDT.gr, drop.self = T, drop.redundant = T, ignore.strand = T)
  out.dt <- data.table(cr1 = tmpDT.gr@elementMetadata$reg[queryHits(xOVR)],
                       cr2 = tmpDT.gr@elementMetadata$reg[subjectHits(xOVR)],
                       int1 = as.integer(substr(x = tmpDT.gr@elementMetadata$reg[queryHits(xOVR)], start = nchar(LIBNAME) + 1, stop = 1000)),
                       int2 = as.integer(substr(x = tmpDT.gr@elementMetadata$reg[subjectHits(xOVR)], start = nchar(LIBNAME) + 1, stop = 1000)))
  out.dt <- rbind(out.dt[int1 <= int2, c("cr1", "cr2"), with = F],
                  data.table(cr1 = out.dt[int1 > int2, cr2], 
                             cr2 = out.dt[int1 > int2, cr1]))
  out.dt[ , reg_con := paste(cr1, cr2, sep = ":")]

  count.dt <- out.dt[ , {read_connection_count = .N; list(read_connection_count = read_connection_count)} , by = "reg_con"]
  
}

count.dt[, c("cr1", "cr2") := tstrsplit(reg_con, ":", fixed = TRUE)]
count.dt <- count.dt[read_connection_count > MINRNUM][order(-read_connection_count)]
count.dt[tmpCOVREG.dt, on = .(cr1 = gene_id), "coordinates_origin" := stringr::str_c(i.seqnames, " ", i.start, " ", i.end, "|", i.strand)]
count.dt[tmpCOVREG.dt, on = .(cr2 = gene_id), "coordinates_target" := stringr::str_c(i.seqnames, " ", i.start, " ", i.end, "|", i.strand)]
count.dt[ , c("cr1", "cr2") := NULL]
count.dt[ , c("sample", "experiment") := .(LIBNAME,  experiment_name)]

# write table
readr::write_csv(x = count.dt, path = file.path(outpath, str_c("smallRNA.connected_regions", experiment_name, LIBNAME, "csv", sep = "."))) 







