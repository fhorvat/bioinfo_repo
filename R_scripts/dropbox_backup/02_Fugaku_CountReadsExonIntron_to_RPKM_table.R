#!/home/students/fhorvat/R/bin/Rscript
### INFO: R Script
### DATE: 05.08.2014
### AUTHOR: Vedran Franke
### PATH: /home/members/vfranke/Projects/PSvoboda_Oocyte_Transcriptome/Scripts/GenomeActivation/ExonAnalysis/CountReadsExonIntron.R
rm(list = ls()); gc()

################################################################################### WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/Fugaku/")

################################################################################### LIBRARIES
lib.path <- "/common/WORK/fhorvat/R_library/vfranke"
source(file.path(lib.path, "FileLoader.R"))
source(file.path(lib.path, "FormatConverters.R"))
source(file.path(lib.path, "BamWorkers.R"))
source(file.path(lib.path, "ScanLib.R"))
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
outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/Fugaku/output/documentation"

inpaths <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/Fugaku/output/Count_ExonIntron_RData"

mm10_gtf_path <- "/common/WORK/fhorvat/reference/mouse/mm10/EnsemblGenes/Ensembl_GRCm38.86.20161128.gtf.gz"

ens_genes_path <- "/common/WORK/fhorvat/reference/mouse/mm10/EnsemblGenes/142007_mm10_EnsemblGenes_Biomart.txt"

tot_counts_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/Fugaku/output/documentation/Fugaku_library_size.txt"

################################################################################### SCRIPT PARAMS
# register workers for parallel computation
registerDoMC(11)

################################################################################### READING TABLES
# read in the repeatmasker
reps <- 
  read_delim("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_repeatMasker_mm10_20170209_all_fields.txt.gz", delim = "\t") %>% 
  dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, repName, repClass, repFamily) %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T)

# read in ENSEMBL annotation with biomart 
ens_annot <- read.table(ens_genes_path, header = TRUE, sep = "\t", stringsAsFactors = F)
ens_annot <- unique(ens_annot[, c("ens.trans.id", "ens.gene.id", "gene.name", "biotype")])

# read in the ensembl annotation
gtf <- read.table(mm10_gtf_path, header=FALSE, sep = "\t", stringsAsFactors = F)
gtf[, 1] <- paste0("chr", gtf[, 1])
gtf <- gtf[!str_detect(gtf[, 1],"NT"),]
gtf[gtf[, 1] == "chrMT", 1] <- "chrM"
gtf <- GffToGRanges(gtf, "exon")

# total counts
tot_counts <- read_delim(file = tot_counts_path, delim = "\t")
tot_counts <- as.integer(tot_counts$library_size) %>% 
  set_names(., tot_counts$ID)

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

# take transcripts with more than one exon
gtf_trans <- gtf_trans[elementNROWS(gtf_trans) > 1]
gtf_trans <- unlist(gtf_trans)
d_val <- data.table(as.data.frame(values(gtf_trans)))
d_val$strand <- as.character(strand(gtf_trans))
d_val[d_val$strand == "+" , `:=`( COUNT = .N , IDX = 1:.N ) , by = transcript_id[strand == "+"]]
d_val[d_val$strand == "-" , `:=`( COUNT = .N , IDX = .N:.1) , by = transcript_id[strand == "-"]]
gtf_trans$ex.num <- d_val$IDX
gtf_trans$ex.tot <- d_val$COUNT
gtf_trans <- split(gtf_trans, as.character(gtf_trans$transcript_id))

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
gtf_range <- unlist(range(gtf_trans))

################################################################################### 
### calculate the intron size
# unlist transcripts to exons, reduce ranges 
red_trans <- reduce(unlist(gtf_trans))

# get exon and total length per gene
fo_int <- data.table(as.matrix(findOverlaps(gtf_range, red_trans)))
fo_int$ex.width <- width(red_trans)[fo_int$subjectHits]
fo_int <- fo_int[,sum(ex.width), by = queryHits]
fo_int$width <- width(gtf_range)[fo_int$queryHits]
setnames(fo_int, 2, "ex.width")

# get length of gene overlaping repeatMasker elements
rep_red <- reduce(reps)
rep_red <- GenomicRanges::setdiff(rep_red, red_trans)
fo_rep_int <- data.table(as.matrix(findOverlaps(gtf_range, rep_red)))
fo_rep_int$rep.width <- width(rep_red)[fo_rep_int$subjectHits]
fo_rep_int <- fo_rep_int[, sum(rep.width), by = queryHits]
setnames(fo_rep_int, 2, "rep.width")

# merge gene/exon/repeatMasker overlaping lengths
int_m <- merge(fo_int, fo_rep_int, by="queryHits", all=TRUE)
int_m$ens.trans.id <- names(gtf_range)[int_m$queryHits]
int_m[is.na(int_m)] <- 0
int_m$queryHits <- NULL

# calculate intron length with/without repeatMasker overlaping lengths
int_m$IntW.Ex <- int_m$width - int_m$ex.width
int_m$IntW.Ex.Rep <- int_m$width - int_m$ex.width - int_m$rep.width
int_m$ex.width <- NULL
int_m$rep.width <- NULL
int_m$width <- NULL

################################################################################### 
# load RData files
rfiles <- unlist(lapply(inpaths, function(x) list.files(x, full.names = T, pattern = "RData")))
rfiles_counts <- rfiles[str_detect(rfiles, "Counts\\.RData")]
rfiles_counts <- rfiles_counts[!str_detect(rfiles_counts, "rep")]
rfiles_counts <- rfiles_counts[!str_detect(rfiles_counts, "PA.SE")]

l_bam <- list()
l_bam <- foreach(i = 1:length(rfiles_counts))%dopar%{
  print(i)
  Assigner(rfiles_counts[i], "l")
  return(l)
}
names(l_bam) <- str_replace(basename(rfiles_counts), ".ExInt_Counts.RData", "")

# read tables from RData files
gene_tabs <- lapply(l_bam, function(x) x[str_detect(names(x), "gene")])
gene_tabs <- lapply(gene_tabs, function(x) x[!str_detect(names(x), "fo_end")])

# change column names 
lapply(names(gene_tabs), function(x)lapply(names(gene_tabs[[x]]), function(y)setnames(gene_tabs[[x]][[y]], 2, paste(x, str_replace(y, "gene.models.", ""), sep = "."))))

# merge tables
gene_tabs <- Reduce(f = function(x, y) merge(x, y, by = "ens.trans.id", all = TRUE), 
                    x = lapply(gene_tabs, function(x) Reduce(f = function(x, y) merge(x, y, by = "ens.trans.id", all = TRUE), x = x)))
gene_tabs[is.na(gene_tabs)] <- 0

################################################################################### 
# create annotation data.frame
gtf_all <- unlist(range(gtf_trans))
sum_width <- sum(width(gtf_trans))
sum_width <- sum_width[sum_width != 0]
ex_num <- elementNROWS(gtf_trans)
ex_num <- ex_num[ex_num != 0]
gene_annot <- data.frame(ens.trans.id = names(gtf_all), 
                         pos = paste(as.character(seqnames(gtf_all)), paste(start(gtf_all), end(gtf_all), sep = "-"), sep = ":"), 
                         ex.width = sum_width, 
                         int.width = width(gtf_all) - sum_width,
                         ex.num = ex_num)
gene_annot <- merge(gene_annot, int_m, by = "ens.trans.id", all = TRUE)

# merges genes annotation with ensembl annotation and loaded genes expression values
m <- merge(gene_annot, ens_annot, by = "ens.trans.id")
gene_d <- merge(m, gene_tabs, by = "ens.trans.id")

# changes colnames
foe_ind <- which(str_detect(colnames(gene_d), "foe"))
setnames(gene_d, foe_ind, str_replace(colnames(gene_d)[foe_ind], "foe", "ex"))
foi_ind <- which(str_detect(colnames(gene_d), "foi"))
setnames(gene_d, foi_ind, str_replace(colnames(gene_d)[foi_ind], "foi", "int"))

# repeats removed
gene_r <- gene_d
gene_r <- cbind(gene_r[, 1:10], gene_r[, str_detect(colnames(gene_r), "(reps.rem.uniq)|(ex.uniq)")])
colnames(gene_r) <- str_replace(colnames(gene_r), ".uniq", "")
colnames(gene_r) <- str_replace(colnames(gene_r), ".rem.uniq", "")
# write.table(gene_r, file.path(outpath, "GeneCount_RepsRemoved_Raw.txt"), row.names=F, col.names=T, quote=F, sep="\t", dec=",")

################################################################################### normalizes gene expression
# gene_m <- gene_d[, !str_detect(colnames(gene_d), "rem")]
# gene_c <- cbind(gene_m[, 1:10], gene_m[, str_detect(colnames(gene_m), "uniq")])
# gene_cnts <- gene_c
# norm_fac <- 1e6 / tot_counts
# 
# # rpm normalization
# samps <- str_replace(names(gene_cnts)[-c(1:10)], "(.int.+)|(.ex.+)", "")
# gene_norm <- cbind(gene_cnts[, 1:10], round(data.frame(t(t(gene_cnts[, -c(1:10)]) * norm_fac[samps])), 1))
# 
# # rpkm normalization
# gene_rpkm <- gene_norm
# ex_ind <- which(str_detect(names(gene_rpkm), "ex.uniq"))
# gene_rpkm[, ex_ind] <- round(gene_rpkm[, ex_ind] * (1e3 / gene_rpkm$ex.width), 2)
# 
# int.ind <- which(str_detect(names(gene_rpkm), "int.uniq"))
# gene_rpkm[, int.ind] <- round(gene_rpkm[, int.ind] * 1e4 / gene_rpkm$int.width, 2)
# gene_rpkm[is.na(gene_rpkm)] <- 0

################################################################################### normalizes for repeats removed
gene_reps <- gene_r
norm_fac <- 1e6 / tot_counts

# rpm normalization
samps <- str_replace(names(gene_reps)[-c(1:10)], "(.int.+)|(.ex.*)", "") 
gene_reps_norm <- cbind(gene_reps[, 1:10], round(data.frame(t(t(gene_reps[, -c(1:10)]) * norm_fac[samps])), 1))

# rpkm normalization
gene_reps_rpkm <- gene_reps_norm
ex_ind <- which(str_detect(names(gene_reps_rpkm), "ex"))[-(1:2)]
gene_reps_rpkm[, ex_ind] <- round(gene_reps_rpkm[, ex_ind] * (1e3 / gene_reps_rpkm$ex.width), 2)
int_ind <- which(str_detect(names(gene_reps_rpkm), "int.reps"))
gene_reps_rpkm[, int_ind] <- round(gene_reps_rpkm[, int_ind] * 1e4 / gene_reps_rpkm$int.width, 2)
gene_reps_rpkm[is.na(gene_reps_rpkm)] <- 0

write.table(gene_reps_rpkm, file.path(outpath, "GeneCount_RepsRemoved_RPKM.txt"), row.names=F, col.names=T, quote=F, sep="\t", dec=",")


