#!/home/students/fhorvat/R/bin/Rscript
### INFO: R Script
### DATE: 4. 4. 2017.  
### AUTHOR: Filip Horvat
### PATH: 
rm(list = ls()); gc()
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

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/CNOT6L_Fugaku_together/documentation"

inpaths <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/CNOT6L_Fugaku_together/RData/01_CountReadsExonIntron"

mm10_gtf_path <- "/common/WORK/fhorvat/reference/mouse/mm10/EnsemblGenes/Ensembl_GRCm38.86.20161128.gtf.gz"

ens_genes_path <- "/common/WORK/fhorvat/reference/mouse/mm10/EnsemblGenes/142007_mm10_EnsemblGenes_Biomart.txt"

tot_counts_path_CNOT6L <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/CNOT6L/documentation/CNOT6L_library_size.txt"
tot_counts_path_Fugaku <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/Fugaku/output/documentation/Fugaku_library_size.txt"

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
registerDoMC(11)

################################################################################### TABLES
################################################################################### 
# read in the repeatmasker
reps <- 
  read_delim("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_repeatMasker_mm10_20170209_all_fields.txt.gz", delim = "\t") %>% 
  dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, repName, repClass, repFamily) %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T)

# read in ENSEMBL annotation with biomart 
ens_annot <- read.table(ens_genes_path, header = TRUE, sep = "\t", stringsAsFactors = F)
ens_annot <- unique(ens_annot[, c("ens.trans.id", "ens.gene.id", "gene.name", "biotype")])

# read in the ensembl annotation
gtf <- read.table(mm10_gtf_path, header = FALSE, sep = "\t", stringsAsFactors = F)
gtf[, 1] <- paste0("chr", gtf[, 1])
gtf <- gtf[!str_detect(gtf[, 1],"NT"),]
gtf[gtf[, 1] == "chrMT", 1] <- "chrM"
gtf <- GffToGRanges(gtf, "exon")

# total counts
sample_table <- rbind(read_delim(file = tot_counts_path_Fugaku, delim = "\t") %>% 
                        dplyr::mutate(name = ID) %>% 
                        dplyr::select(ID, name, library_size), 
                      read_delim(file = tot_counts_path_CNOT6L, delim = "\t"))

tot_counts <- as.integer(sample_table$library_size) %>% 
  set_names(., sample_table$name)

norm_counts <- 1e6 / tot_counts

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

################################################################################### 
# load RData files as data.frame
rfiles <- unlist(lapply(inpaths, function(x) list.files(x, full.names = T, pattern = "RData")))

# get file names
rfiles_names <-
  tibble(rfiles = basename(rfiles)) %>%
  dplyr::mutate(ID = ifelse(str_detect(rfiles, "s_"),
                            str_replace(rfiles, ".ExInt_Counts.RData", ""), 
                            str_replace(rfiles, "_.*", ""))) %>% 
  left_join(., sample_table, by = "ID") %$% 
  name

# load files in parallel 
gene_tabs <- list()
gene_tabs <- foreach(i = 1:length(rfiles))%dopar%{
  print(i)
  Assigner(rfiles[i], "l")
  return(l)
}

# set names
names(gene_tabs) <- rfiles_names
  
# change column names 
lapply(names(gene_tabs), function(x) lapply(names(gene_tabs[[x]]), function(y) setnames(gene_tabs[[x]][[y]], 2, paste(x, str_replace(y, "gene.models.", ""), sep = "."))))

# merge tables
gene_tabs <- Reduce(f = function(x, y) merge(x, y, by = "ens.trans.id", all = TRUE), 
                    x = lapply(gene_tabs, function(x) Reduce(f = function(x, y) merge(x, y, by = "ens.trans.id", all = TRUE), x = x)))
gene_tabs[is.na(gene_tabs)] <- 0

################################################################################### 
# merges genes annotation with ensembl annotation and loaded genes expression values
m <- merge(gene_annot, ens_annot, by = "ens.trans.id")
gene_reps <- merge(m, gene_tabs, by = "ens.trans.id")

# changes colnames
foe_ind <- which(str_detect(colnames(gene_reps), "foe"))
setnames(gene_reps, foe_ind, str_replace(colnames(gene_reps)[foe_ind], "foe", "ex"))
foi_ind <- which(str_detect(colnames(gene_reps), "foi"))
setnames(gene_reps, foi_ind, str_replace(colnames(gene_reps)[foi_ind], "foi", "int"))

# repeats removed
gene_reps <- gene_reps
gene_reps <- cbind(gene_reps[, 1:10], gene_reps[, str_detect(colnames(gene_reps), "(reps_rem_uniq)|(ex_uniq)")])
colnames(gene_reps) <- str_replace(colnames(gene_reps), "_uniq", "")
colnames(gene_reps) <- str_replace(colnames(gene_reps), "_rem_uniq", "")

################################################################################### normalizes for repeats removed
# rpm normalization
samps <- str_replace(names(gene_reps)[-c(1:10)], "(.int.+)|(.ex.*)", "") 
gene_reps_rpkm <- cbind(gene_reps[, 1:10], round(data.frame(t(t(gene_reps[, -c(1:10)]) * norm_counts[samps])), 1))

# rpkm normalization
ex_ind <- which(str_detect(names(gene_reps_rpkm), "ex"))[-(1:2)]
gene_reps_rpkm[, ex_ind] <- round(gene_reps_rpkm[, ex_ind] * (1e3 / gene_reps_rpkm$ex.width), 2)
int_ind <- which(str_detect(names(gene_reps_rpkm), "int_reps"))
gene_reps_rpkm[, int_ind] <- round(gene_reps_rpkm[, int_ind] * 1e4 / gene_reps_rpkm$int.width, 2)
gene_reps_rpkm[is.na(gene_reps_rpkm)] <- 0

# write table
write.table(gene_reps_rpkm, file.path(outpath, "GeneCount_RepsRemoved_RPKM.txt"), row.names = F, col.names = T, quote = F, sep = "\t", dec = ",")

################################################################################### average of RPKM in samples
# get sample table
rpkm_table <- gene_reps_rpkm[, -c(1:10)]

# get indices of samples which are in replicates
rep_ind <- which(str_detect(colnames(rpkm_table), "X"))

# get indices of samples which are not in replicates
no_rep_ind <- setdiff(1:ncol(rpkm_table), rep_ind)

# get unique sample names
sample_names_unique <- 
  str_replace(colnames(rpkm_table)[rep_ind], ".int_reps_rem|.ex", "") %>% 
  str_replace(., "_X.*", "") %>% 
  unique()

gene_reps_rpkm_average <- 
  lapply(sample_names_unique, function(x){
    
    # get one sample columns
    rpkm_table_average <- rpkm_table[, str_detect(colnames(rpkm_table), x)]

    # order table (exon - intron)
    rpkm_table_average <- rpkm_table_average[, c(which(str_detect(colnames(rpkm_table_average), "ex")), 
                                                 which(str_detect(colnames(rpkm_table_average), "int")))]
    
    # calculate average RPKM expression for sample (exon / intron)
    rpkm_table_average <- 
      cbind(apply(X = rpkm_table_average[, 1:(ncol(rpkm_table_average) / 2)], MARGIN = 1, FUN = mean), 
            apply(X = rpkm_table_average[, ((ncol(rpkm_table_average) / 2) + 1):ncol(rpkm_table_average)], MARGIN = 1, FUN = mean)) %>% 
      as.data.frame(.) %>% 
      round(., 2) %>% 
      set_colnames(c(str_c(x, ".ex", sep = ""), 
                     str_c(x, ".int_reps_rem", sep = "")))
    
    return(rpkm_table_average)
    
  }) %>% 
  bind_cols(.)

# bind with rest of data
gene_reps_rpkm_average <- 
  cbind(gene_reps_rpkm[, c(1:10)], 
        rpkm_table[no_rep_ind], 
        gene_reps_rpkm_average)

# write as table
write.table(gene_reps_rpkm_average, file.path(outpath, "GeneCount_RepsRemoved_RPKM_average.txt"), row.names = F, col.names = T, quote = F, sep = "\t", dec = ",")
