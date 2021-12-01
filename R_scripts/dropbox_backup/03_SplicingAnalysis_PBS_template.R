#!/home/students/fhorvat/R/bin/Rscript
### INFO: R Script
### DATE: 5. 4. 2017.  
### AUTHOR: Filip Horvat
### PATH: 
rm(list = ls()); gc()
# options(bitmapType = 'cairo')

################################################################################### WORKING DIRECTORY
################################################################################### 
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads")

################################################################################### LIBRARIES
################################################################################### 
library(doMC)
library(GenomicAlignments)
library(genomation)
library(data.table)
library(stringr)
library(readr)
library(dplyr)
library(magrittr)

################################################################################### PATH VARIABLES
################################################################################### 
lib_path <- "/common/WORK/fhorvat/R_library/vfranke"

r_outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/CNOT6L_Fugaku_together/RData/03_SplicingAnalysis"
dir.create(r_outpath, showWarnings = FALSE)

mm10_gtf_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl_GRCm38.86.20161128.gtf.gz"

genes_sel_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/CNOT6L_Fugaku_together/documentation/GeneCount_RepsRemoved_RPKM_average.txt"

# file <- %FILE

################################################################################### SOURCE FILES
################################################################################### 
source(file.path(lib_path, "FileLoader.R"))
source(file.path(lib_path, "FormatConverters.R"))
source(file.path(lib_path, "BamWorkers.R"))
source(file.path(lib_path, "ScanLib.R"))

################################################################################### FUNCTIONS
################################################################################### 

################################################################################### SCRIPT PARAMS
################################################################################### 
# register workers for parallel computation
registerDoMC(21)

################################################################################### TABLES
################################################################################### 
# repeatMasker table
reps <- 
  read_delim("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_repeatMasker_mm10_20170209_all_fields.txt.gz", delim = "\t") %>% 
  dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, repName, repClass, repFamily) %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T)

# reads in the ensembl annotation
gtf <- read.table(mm10_gtf_path, header = FALSE, sep = "\t", stringsAsFactors = F)
gtf[, 1] <- paste0("chr", gtf[, 1])
gtf <- gtf[!str_detect(gtf[, 1],"NT"),]
gtf[gtf[, 1] == "chrMT", 1] <-"chrM"
gtf <- GffToGRanges(gtf, "exon")

################################################################################### MAIN CODE
################################################################################### 
# get unique gene_id-transcript_id combinations
gids <- unique(values(gtf)[c("gene_id", "transcript_id")])

# create ENSEMBL transcripts ranges, takes one transcript for each gene (the one with the most exons)
gtf_trans <- split(gtf, gtf$transcript_id)
gtf_trans <- gtf_trans[order(-elementNROWS(gtf_trans))]
gtf_trans <- gtf_trans[!duplicated(gids$gene_id[match(names(gtf_trans), gids$transcript_id)])]

# take transcripts with single exon
gtf_trans_single <- gtf_trans[elementNROWS(gtf_trans) == 1]
gtf_trans_single <- unlist(gtf_trans_single)
values(gtf_trans_single) <- gtf_trans_single$transcript_id
names(mcols(gtf_trans_single)) <- "transcript_id"
gtf_trans_single$ex.num <- 1
gtf_trans_single$ex.tot <- 1
gtf_trans_single <- split(gtf_trans_single, gtf_trans_single$transcript_id)

# take transcripts with more than one exon
gtf_trans <- gtf_trans[elementNROWS(gtf_trans) > 1]
gtf_trans <- reduce(gtf_trans, min.gapwidth = 100)
gtf_trans <- unlist(gtf_trans)
gtf_trans$transcript_id <- names(gtf_trans)
d_val <- data.table(as.data.frame(values(gtf_trans)))
d_val$strand <- as.character(strand(gtf_trans))
d_val[d_val$strand == "+" , `:=`( COUNT = .N , IDX = 1:.N ), by = transcript_id[strand == "+"]]
d_val[d_val$strand == "-" , `:=`( COUNT = .N , IDX = .N:.1), by = transcript_id[strand == "-"]]
gtf_trans$ex.num <- d_val$IDX
gtf_trans$ex.tot <- d_val$COUNT
gtf_trans <- split(gtf_trans, as.character(gtf_trans$transcript_id))

# get exons
gtf_ex <- unlist(gtf_trans)

################################################################################### 
# reads in the selected protein coding genes
genes_sel <- 
  read_delim(file = genes_sel_path, delim = "\t", locale = locale(decimal_mark = ",")) %>% 
  dplyr::filter(biotype == "protein_coding") %>% 
  as.data.frame(.)

# get genes which are not expressed in exons/introns in MII stage and which are expressed in 1-cell stage introns
genes_sel_mii <- 
  genes_sel %>% 
  dplyr::filter(s_MII_WT.int_reps_rem == 0, s_MII_WT.ex == 0, s_1C_WT.int_reps_rem != 0)

# create list of all those subsets
l_genes <- list(genes_sel_mii = genes_sel_mii)

################################################################################### 
# constructs the splicing coordinates
# splice donor
splice <- unlist(gtf_trans)
splice_don <- GenomicRanges::resize(splice, fix = "end", width = 1)
splice_don <- splice_don[splice_don$ex.num != splice_don$ex.tot]
splice_don <- resize(splice_don, width = 10, fix = "start")
splice_don <- unique(splice_don)
splice_don <- splice_don[countOverlaps(splice_don, splice_don) == 1]

# splice acceptor
splice_acc <- GenomicRanges::resize(splice, fix = "start", width = 1)
splice_acc <- splice_acc[splice_acc$ex.num != 1]	
splice_acc <- resize(splice_acc, width = 10, fix = "end")
splice_acc <- unique(splice_acc)
splice_acc <- splice_acc[countOverlaps(splice_acc, splice_acc) == 1]

################################################################################### 
# unite transcripts with one and more exons
gtf_genes <- c(gtf_trans, gtf_trans_single)

# get full range of those transcripts
gtf_range <- unlist(range(gtf_genes))

# get introns 
gtf_int <- GenomicRanges::setdiff(gtf_range, gtf_ex)

################################################################################### 
# loops through the file and calculates the statistics
name <- BamName(file)
print(name)
chrs <- chrFinder(file)
chrs <- chrs[chrs$chr != "chrM",]

################################################################################### 
# loops through the file and calculates the statistics
name <- BamName(file)
print(name)
chrs <- chrFinder(file)
chrs <- chrs[chrs$chr != "chrM",]

# for each chromosome in parallel 
reads_stat <- foreach(chr = chrs$chr) %dopar% {
  
  # print chromosome name
  print(chr)
  
  ################################################################################### 
  # read in the reads
  cat("Reading files...\n")
  which_ranges <- GRanges(chr, IRanges(1, chrs$chrlen[chrs$chr == chr]))
  if(str_detect(name, "SE")){
    bam <- readGAlignments(file, param = ScanBamParam(which = which_ranges, tag = "NH"), use.names = TRUE)
    bams <- bam
  }else{
    bam <- suppressWarnings(readGAlignmentPairs(file, param = ScanBamParam(which = which_ranges), use.names = TRUE))
    bams <- readGAlignments(file, param = ScanBamParam(which = which_ranges, tag = "NH"), use.names = TRUE)
  }
  
  ################################################################################### 
  # create grangesList for reads, get reads which are uniquely mapped
  gbam <- grglist(bam)
  tab <- unique(data.table(names = names(bams), V1 = values(bams)$NH))
  tab <- tab[tab$names %in% names(bam)]
  tab <- tab[tab$V1 == 1]
  rm(bam, bams); gc()
  
  ################################################################################### 
  # count overlaps between reads and exons/introns
  cat("Selecting reads...\n")
  co_ex <- countOverlaps(gbam, gtf_ex, ignore.strand = TRUE)
  co_in <- countOverlaps(gbam, gtf_int, ignore.strand = TRUE, minoverlap = 15)
  
  # get indices of reads which overlap more than one exon or at least one exon and intron
  co_ind <- which(co_ex > 1 | (co_ex >= 1 & co_in >= 1))
  
  ################################################################################### 
  # table of reads
  ######### map to transcripts
  # find overlaps between reads and full transcript ranges
  fot <- data.table(as.matrix(findOverlaps(gbam, gtf_range, ignore.strand = TRUE)))
  fot$id <- names(gbam)[fot$queryHits]
  
  ######### remove reads overlaping more transcripts
  fon <- fot[, length(unique(subjectHits)), by = id]
  fon <- fon[fon$V1 == 1]
  fot <- fot[fot$id %in% fon$id]
  
  ######### transcript ID
  # to each read add transcript ID to which that read is mapped to
  fot$ens.trans.id <- names(gtf_range)[fot$subjectHits]
  
  ######### expression
  # for each transcript label if it belongs to :
  #  - genes which are not expressed in exons/introns in MII stage and which are expressed in introns in 1-cell stage 
  fot$genes_sel_mii <- fot$ens.trans.id %in% l_genes[["genes_sel_mii"]]$ens.trans.id
  
  ######### distance
  # add distance between two reads in one paired-end fragment  
  fot$dist <- max(start(gbam)[fot$queryHits] - min(end(gbam)[fot$queryHits]))
  # add length of the fragment
  fot$dist.2 <- max(end(gbam)[fot$queryHits]) - min(start(gbam)[fot$queryHits])
  
  ######### repeatMasker
  # label if read is overlapping and range from repeatMasker
  forep <- data.table(as.matrix(findOverlaps(gbam, reps, ignore.strand = TRUE)))
  fot$reps <- fot$queryHits %in% forep$queryHits
  
  ######### filter reads
  # take reads:
  #  - mapped to more than one exon or at least one exon and intron in this .bam file
  #  - uniquely mapped
  #  - not mappped to repeat from repeatMasker
  #  - unique reads
  fot <- fot[fot$queryHits %in% co_ind]
  fot <- fot[fot$id %in% tab$names]
  fot <- fot[!fot$reps]
  fot$ens.trans.id <- NULL
  fot$subjectHits <- NULL
  fot <- unique(fot)
  
  ######### exon/intron 
  # label if read overlaps more than one exon
  fot$ex <- fot$queryHits %in% which(co_ex > 1)
  # label if read overlaps one or more introns
  fot$int <- fot$queryHits %in% which(co_in > 0)
  # if read overlaps any introns set exon-overlap label to FALSE
  fot$ex[fot$int > 0] <- FALSE
  
  ######### return
  # return as list 
  return(fot)
  
} %>% 
  rbindlist(.)

# save data
cat("Saving data...\n")
save(reads_stat, file = file.path(r_outpath, paste(name, "SpliceDonAcc.RData", sep = ".")))






