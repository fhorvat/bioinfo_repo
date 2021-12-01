#!/home/students/fhorvat/R/bin/Rscript
### INFO: R Script
### DATE: 05.08.2014
### AUTHOR: Vedran Franke
### PATH: /home/members/vfranke/Projects/PSvoboda_Oocyte_Transcriptome/Scripts/GenomeActivation/ExonAnalysis/PBS_CountReadsExonIntron.Temp.R
rm(list = ls()); gc()

################################################################################### WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/Fugaku/")

################################################################################### LIBRARIES
lib_path <- "/common/WORK/fhorvat/R_library/vfranke"
source(file.path(lib_path, "FileLoader.R"))
source(file.path(lib_path, "FormatConverters.R"))
source(file.path(lib_path, "BamWorkers.R"))
source(file.path(lib_path, "ScanLib.R"))
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
# file <- %FILE
file <- '/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_1cell.WE/s_1cell.WE.bam'

r_outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/Fugaku/output/Count_ExonIntron_RData"
# dir.create(r_outpath, showWarnings = FALSE)

mm10_gtf_path <- "/common/WORK/fhorvat/reference/mouse/mm10/EnsemblGenes/Ensembl_GRCm38.86.20161128.gtf.gz"

ooreg_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/Fugaku/output/documentation/OocyteRegions.lower.5.bed"

################################################################################### SCRIPT PARAMS
# register workers for parallel computation
registerDoMC(21)

################################################################################### READING TABLES
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

# read in oocyte regions
oo_regs <- readGeneric(ooreg_path)

################################################################################### MAIN CODE
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

###################################################################################
### get ranges of areas 10000 bases downstream of transcripts in 1kb bins
# get end of gene
gtf_end <- resize(gtf_range, width = 1, fix = "end")
gtf_end <- resize(gtf_end, width = 1000, fix = "start")
gtf_end$ens.trans.id <- names(gtf_end)

# add 1kb bins downstream of gene end
l_end <- list()
l_end[[1]] <- gtf_end
for(i in 2:10){
  l_end[[i]] <- flank(l_end[[i - 1]], width = 1000, start = FALSE)
}

# get GRanges, order by start
gend <- unlist(GRangesList(l_end))
gend <- gend[order(start(gend))]

# create data.table, add strand data
g_val <- data.table(as.data.frame(values(gend)))
g_val$strand <- as.character(strand(gend))

# enumerate 1kb bins downstream of gene
g_val[g_val$strand == "+" , `:=`( COUNT = .N , IDX = 1:.N ) , by = ens.trans.id[strand == "+"]]
g_val[g_val$strand == "-" , `:=`( COUNT = .N , IDX = .N:.1) , by = ens.trans.id[strand == "-"]]

# add index of 1kb bin as new column in data.table
gend$reg.num <- g_val$IDX

###################################################################################
### construct the splicing coordinates
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
# unite transcripts with one and more than one exons
gtf_trans <- c(gtf_trans, gtf_trans_single)

# get ranges of whole gene for that transcript (exons + introns)
# gtf_range <- unlist(range(gtf_trans))

###################################################################################
# loops through all files

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
  
  # add weight based on number of positions in genome read is mapped to
  foe$weight <- round(1 / tab$V1[match(foe$id, tab$names)], 2)
  
  # add nweight based on number of transcripts read is mapped to
  foe$nweight <- round(1 / foe$N, 2)
  
  # for each transcript sum how many uniquely mapped reads are mapped to that transcript
  foe_uniq <- foe[foe$id %in% tab$names[tab$V1 == 1] & foe$N == 1, length(id), by = "ens.trans.id"]
  
  # for each transcript sum weight and nweight values of reads mapped to that transcript
  foe_weig <- foe[, sum(weight * nweight), by = "ens.trans.id"]
  ###
  
  ### counts reads mapped to introns (ones which are not mapped to any exons)
  fogr <- data.table(as.matrix(findOverlaps(gbam, gtf_range, ignore.strand = TRUE)))
  foi <- fogr[!fogr$queryHits %in% foe$queryHits]
  foi$ens.trans.id <- names(gtf_trans)[foi$subjectHits]
  foi_n <- foi[, .N, by = queryHits]
  foi <- merge(foi, foi_n, by = "queryHits", all = TRUE)
  foi$id <- names(gbam)[foi$queryHits]
  foi$weight <- round(1 / tab$V1[match(foi$id, tab$names)], 2)
  foi$nweight <- round(1 / foi$N, 2)
  foi_uniq <- foi[foi$id %in% tab$names[tab$V1 == 1] & foi$N == 1, length(id), by="ens.trans.id"]
  foi_weig <- foi[, sum(weight * nweight), by = "ens.trans.id"]
  ###
  
  ###################################################################################
  cat("Counting exon donors and acceptors...\n")
  ### find the reads mapped to splice junctions
  # splice donors
  fo_don <- data.table(as.matrix(findOverlaps(splice_don, gbam, ignore.strand = TRUE, type = "within")))
  fo_don$id <- names(gbam)[fo_don$subjectHits]
  
  # splice acceptors
  fo_acc <- data.table(as.matrix(findOverlaps(splice_acc, gbam, ignore.strand = TRUE, type = "within")))
  fo_acc$id <- names(gbam)[fo_acc$subjectHits]
  ###
  
  ###################################################################################
  cat("Counting repeats...\n")
  ### find the reads mapped to repeats from repeatMasker table
  fo_rep <- data.table(as.matrix(findOverlaps(reps, gbam, ignore.strand=TRUE)))
  fo_rep$id <- names(gbam)[fo_rep$subjectHits]
  
  # add weight based on number of positions in genome read is mapped to
  fo_rep$weight <- round(1 / tab$V1[match(fo_rep$id, tab$names)], 2)
  
  # for each repeat sum how many uniquely mapped reads are mapped to that repeat
  co_reps <- fo_rep[, length(subjectHits), by = queryHits]
  
  # for each repeat sum weights of reads mapped to that reapeat
  co_weig <- fo_rep[, sum(weight), by = queryHits]
  ###
  
  ###################################################################################
  ### find oocyte reads (reads which are mapped to regions which also have mapped reads in MII stage)
  foe_oo <- data.table(as.matrix(findOverlaps(gbam, oo_regs, ignore.strand = TRUE)))
  foe_oo$names <- names(gbam)[foe_oo$queryHits]
  
  ################################################################################### remove reads
  ### remove oocyte reads from reads mapped to introns
  foi_oo_rem <- foi[!foi$id %in% foe_oo$names]
  
  # for each transcript sum how many uniquely mapped reads not mapped to oocyte regions are mapped to that transcript
  foi_oo_rem_uniq <- foi_oo_rem[foi_oo_rem$id %in% tab$names[tab$V1 == 1] & foi_oo_rem$N == 1, length(id), by = "ens.trans.id"]
  
  # for each transcript create sum weight and nweight of reads not mapped to oocyte regions
  foi_oo_rem_weig <- foi_oo_rem[, sum(weight * nweight), by = "ens.trans.id"]
  ###
  
  ### remove reads mapped to repeatMasker regions from reads mapped to introns
  foi_reps_rem <- foi[!foi$id %in% fo_rep$id]
  
  # for each transcript sum how many uniquely mapped reads not mapped to repeatMasker regions are mapped to that transcript
  foi_reps_rem_uniq <- foi_reps_rem[foi_reps_rem$id %in% tab$names[tab$V1 == 1] & foi_reps_rem$N == 1, length(id), by = "ens.trans.id"]
  
  # for each transcript sum weight and nweight values of reads not mapped to repeatMasker regions
  foi_reps_rem_weig <- foi_reps_rem[, sum(weight * nweight), by = "ens.trans.id"]
  ###
  
  ###################################################################################
  ### find reads mapped to 1kb bins 10kb downstream of gene
  fo_end <-  data.table(as.matrix(findOverlaps(gbam, gend, ignore.strand = TRUE)))
  fo_end$id <- names(gbam)[fo_end$queryHits]
  
  # remove reads mapped to repeatMasker regions
  fo_end <- fo_end[!fo_end$id %in% fo_rep$id]
  
  # add ID of transcripts from which are reads mapped downstream
  fo_end$ens.trans.id <- gend$ens.trans.id[fo_end$subjectHits]
  
  # add index of bin downstream from gene to which the read is mapped to
  fo_end$reg.num <- gend$reg.num[fo_end$subjectHits]
  
  # order reads based on index of bin reads are mapped to
  fo_end <- fo_end[order(fo_end$reg.num)]
  
  # take only unique combinations of transcripts and read ID mapped to them
  fo_end <- fo_end[!duplicated(paste(fo_end$ens.trans.id, fo_end$id))]
  
  # take only reads which are mapped to unique position in genome
  fo_end <- fo_end[fo_end$id %in% tab$names[tab$V1 == 1]]
  
  # count how many reads are mapped to each 1kb bin in each transcript in this sample
  fo_end <- fo_end[, length(queryHits), by = c("ens.trans.id", "reg.num")]
  setnames(fo_end, 3, name)
  
  ###################################################################################
  ### table for read statistics
  # is read mapped to exon?
  tab$exon <- tab$names %in% foe$id
  
  # is read mapped to intron?
  tab$intr <- tab$names %in% foi$id
  
  # is read mapped to repeatMasker region?
  tab$reps <- tab$names %in% fo_rep$id
  
  # is read mapped to splice donor?
  tab$sp.don <- tab$names %in% fo_don$id
  
  # is read mapped to splice acceptor?
  tab$sp.acc <- tab$names %in% fo_acc$id
  ###
  
  ###################################################################################
  # return list which consists of:
  # - gene models tables = read counts/weights of reads:
  #   - uniquely mapped to exons/introns (foe_uniq/weig, foi_uniq/weig)
  #   - uniquely mapped to introns of non-oocyte regions (foi_oo_rem_uniq/weig)
  #   - uniquely mapped to introns of non-repeatMasker regoins (foi_reps_rem_uniq/weig)
  #   - uniquely mapped to 1kb bins 10kb from gene ends (fo_end)
  # - repeats tables = read counts/weights of reads:
  #   - uniqeuly mapped to repeatMasker regions (co_reps/weig)
  # - tabs = table of read statistics, is each read mapped to (tabs):
  #   - exon
  #   - intron
  #   - repeatMasker region
  #   - splice donor
  #   - splice acceptor
  
  return(list(gene.models = list(tables = list(foe_uniq = foe_uniq,
                                               foi_uniq = foi_uniq,
                                               foe_weig = foe_weig,
                                               foi_weig = foi_weig,
                                               foi_oo_rem_uniq = foi_oo_rem_uniq,
                                               foi_oo_rem_weig = foi_oo_rem_weig,
                                               foi_reps_rem_uniq = foi_reps_rem_uniq,
                                               foi_reps_rem_weig = foi_reps_rem_weig,
                                               fo_end = fo_end)),
              reps = list(tables = list(co_reps = co_reps,
                                        co_weig = co_weig)),
              tabs = list(reads = tab)))
  
}

###################################################################################
cat("Saving data...\n")
# extract counts/weights of reads mapped to gene models and repeat tables for each chromosome
l_table <- lapply(cnts, function(x) unlist(lapply(x, function(y) y [["tables"]]), recursive = FALSE))

# bind together each list into data.frames, save as RData
d_table <- lapply(names(l_table[[1]]), function(x) rbindlist(lapply(l_table, "[[", x)))
names(d_table) <- names(l_table[[1]])
save(d_table, file = file.path(r_outpath, paste(name, "ExInt_Counts.RData", sep = ".")))

# extract and save read statistics
d_reads <- rbindlist(lapply(cnts, function(x) x$tabs$reads))
save(d_reads, file = file.path(r_outpath, paste(name, "ExInt_Table.RData", sep = ".")))

