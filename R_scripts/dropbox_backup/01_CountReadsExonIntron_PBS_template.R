#!/home/students/fhorvat/R/bin/Rscript
### INFO: R Script
### DATE: 4. 4. 2017. 
### AUTHOR: Filip Horvat
### PATH: 
rm(list = ls()[!str_detect(ls(), "^gtf$|^reps$")]); gc()
# options(bitmapType = 'cairo')

################################################################################### WORKING DIRECTORY
################################################################################### 
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads")

################################################################################### LIBRARIES
################################################################################### 
library(data.table)
library(stringr)
library(doMC)
library(GenomicAlignments)
library(genomation)
library(reshape2)
library(readr)
library(dplyr)
library(magrittr)

################################################################################### PATH VARIABLES
################################################################################### 
lib_path <- "/common/WORK/fhorvat/R_library/vfranke"

r_outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/CNOT6L_Fugaku_together/RData/01_CountReadsExonIntron"
dir.create(r_outpath, showWarnings = FALSE)

mm10_gtf_path <- "/common/WORK/fhorvat/reference/mouse/mm10/EnsemblGenes/Ensembl_GRCm38.86.20161128.gtf.gz"

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
# read in the repeatmasker
reps <-
  read_delim("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_repeatMasker_mm10_20170209_all_fields.txt.gz", delim = "\t") %>%
  dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, repName, repClass, repFamily) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = T)

# read in the ensembl annotation
gtf <- read.table(mm10_gtf_path, header=FALSE, sep = "\t", stringsAsFactors = F)
gtf[, 1] <- paste0("chr", gtf[, 1])
gtf <- gtf[!str_detect(gtf[, 1],"NT"),]
gtf[gtf[, 1] == "chrMT", 1] <- "chrM"
gtf <- GffToGRanges(gtf, "exon")

################################################################################### MAIN CODE
################################################################################### 
# create ENSEMBL transcripts ranges, takes one transcript for each gene (the one with the most exons)
gids <- unique(values(gtf)[c("gene_id", "transcript_id")])
gtf_trans <- split(gtf, gtf$transcript_id)
gtf_trans <- gtf_trans[order(-elementNROWS(gtf_trans))]
gtf_trans <- gtf_trans[!duplicated(gids$gene_id[match(names(gtf_trans), gids$transcript_id)])]

# take transcripts with single exon
gtf_trans_single <- gtf_trans[elementNROWS(gtf_trans) == 1]
gtf_trans_single <- unlist(gtf_trans_single)
gtf_trans_single$ex.num <- 1
gtf_trans_single$ex.tot <- 1
gtf_trans_single <- split(gtf_trans_single, gtf_trans_single$transcript_id)

# take transcripts with more than one exon, count and enumerate exons in each transcript
gtf_trans <- gtf_trans[elementNROWS(gtf_trans) > 1]
gtf_trans <- unlist(gtf_trans)
d_val <- data.table(as.data.frame(values(gtf_trans)))
d_val$strand <- as.character(strand(gtf_trans))
d_val[d_val$strand == "+" , `:=`( COUNT = .N , IDX = 1:.N ) , by = transcript_id[strand == "+"]]
d_val[d_val$strand == "-" , `:=`( COUNT = .N , IDX = .N:.1) , by = transcript_id[strand == "-"]]
gtf_trans$ex.num <- d_val$IDX
gtf_trans$ex.tot <- d_val$COUNT
gtf_trans <- split(gtf_trans, as.character(gtf_trans$transcript_id))

# get range of transcripts (= range of gene which is transcribed into that transcript)
gtf_range <- unlist(range(gtf_trans))

# unite transcripts with one and more than one exons
gtf_trans <- c(gtf_trans, gtf_trans_single)

###################################################################################
# get file name and chromosomes in .bam file
name <- BamName(file)
print(name)
chrs <- chrFinder(file)
chrs <- chrs[chrs$chr != "chrM",]

# loops through the file in parallel and calculates the statistics
cnts <- foreach(chr = chrs$chr)%dopar%{
  
  print(chr)
  
  # reads in the reads
  cat("Reading files...\n")
  which_ranges <- GRanges(chr, IRanges(1, chrs$chrlen[chrs$chr == chr]))
  
  if(str_detect(name, "SE")){
    bam <- readGAlignments(file, param = ScanBamParam(which = which_ranges, tag = "NH"), use.names = TRUE)
    bams <- bam
  }else{
    bam <- suppressWarnings(readGAlignmentPairs(file, param = ScanBamParam(which = which_ranges), use.names = TRUE))
    bams <- readGAlignments(file, param = ScanBamParam(which = which_ranges, tag = "NH"), use.names = TRUE)
  }
  
  # make GRanges List for reads
  gbam <- grglist(bam)
  tab <- unique(data.table(names = names(bams), V1 = values(bams)$NH))
  tab <- tab[tab$names %in% names(bam)]
  rm(bam, bams); gc()
  
  ###################################################################################
  cat('Counting reads in models...\n')
  
  ###  count reads mapped to exons
  # find overlaps between reads and transcripts
  foe <- data.table(as.matrix(findOverlaps(gbam, gtf_trans, ignore.strand = TRUE)))
  
  # for each read add ID of transcript read is mapped to
  foe$ens.trans.id <- names(gtf_trans)[foe$subjectHits]
  
  # for each read count number of transcripts read is mapped to (almost always = 1)
  foe_n <- foe[, .N, by = queryHits]
  foe <- merge(foe, foe_n, by = "queryHits", all = TRUE)
  
  # add ID of each read
  foe$id <- names(gbam)[foe$queryHits]
  
  # for each transcript sum how many uniquely mapped reads are mapped to that transcript
  foe_uniq <- foe[foe$id %in% tab$names[tab$V1 == 1] & foe$N == 1, length(id), by = "ens.trans.id"]
  ###
  
  ### counts reads mapped to introns (ones which are not mapped to any exons)
  fogr <- data.table(as.matrix(findOverlaps(gbam, gtf_range, ignore.strand = TRUE)))
  foi <- fogr[!fogr$queryHits %in% foe$queryHits]
  foi$ens.trans.id <- names(gtf_trans)[foi$subjectHits]
  foi_n <- foi[, .N, by = queryHits]
  foi <- merge(foi, foi_n, by = "queryHits", all = TRUE)
  foi$id <- names(gbam)[foi$queryHits]
  foi_uniq <- foi[foi$id %in% tab$names[tab$V1 == 1] & foi$N == 1, length(id), by="ens.trans.id"]
  ###
  
  ### find the reads mapped to repeats from repeatMasker table
  fo_rep <- data.table(as.matrix(findOverlaps(reps, gbam, ignore.strand = TRUE)))
  fo_rep$id <- names(gbam)[fo_rep$subjectHits]
  
  # remove reads mapped to repeatMasker regions from reads mapped to introns
  foi_reps_rem <- foi[!foi$id %in% fo_rep$id]
  
  # for each transcript sum how many uniquely mapped reads not mapped to repeatMasker regions are mapped to that transcript
  foi_reps_rem_uniq <- foi_reps_rem[foi_reps_rem$id %in% tab$names[tab$V1 == 1] & foi_reps_rem$N == 1, length(id), by = "ens.trans.id"]
  ###
  
  ###################################################################################
  # return list which consists of:
  # - gene models tables = counts of reads:
  #   - uniquely mapped to exons (foe_uniq)
  #   - uniquely mapped to introns (foi_uniq)
  #   - uniquely mapped to introns of non-repeatMasker regions (foi_reps_rem_uniq)
  
  return(gene.models = list(foe_uniq = foe_uniq,
                            foi_uniq = foi_uniq, 
                            foi_reps_rem_uniq = foi_reps_rem_uniq))
  
}

###################################################################################
# bind together each list into data.frames
cat("Saving data...\n")
d_table <- 
  lapply(names(cnts[[1]]), function(x) rbindlist(lapply(cnts, "[[", x))) %>% 
  set_names(., names(cnts[[1]]))

# save as .RData
save(d_table, file = file.path(r_outpath, paste(name, "ExInt_Counts.RData", sep = ".")))

