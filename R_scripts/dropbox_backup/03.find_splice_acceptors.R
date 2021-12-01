### INFO: 
### DATE: Mon Mar 19 14:39:40 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/Ago2_repeats_evolution/splice_sites")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicFeatures)
library(GenomicRanges)
library(Rsamtools)
library(Biostrings)
library(RWebLogo)
library(seqLogo)
library(ggbio)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(BSgenome.Mmusculus.UCSC.mm10)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "gtfToGRanges.R"))
source(file.path(lib_path, "mutate_cond.R"))
wideScreen()

######################################################## FUNCTIONS
# calculate proportion
proportion <- function(x){
  rs <- sum(x)
  return(x / rs)
}

######################################################## PATH VARIABLES
# set outpath
outpath <- getwd()

# gtf path
gtf_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/ensembl.Mus_musculus.GRCm38.89.20180305.UCSCnames.clean.gtf.gz"

# gtf info
gtf_info_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/ensembl.Mus_musculus.GRCm38.89.20180305.UCSCnames.clean.geneInfo.csv"

# mm10 fasta path 
fasta_path <- "/common/WORK/fhorvat/reference/mouse/mm10/UCSC/mm10.fa"

# sequences path
seq_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/Ago2_repeats_evolution/sequences/Ago2.intron3_4.mammals.sense.fasta"

# bam path
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/Ago2_repeats_evolution/expression/s_GV.WE.Ago2.bam"

######################################################## READ DATA
# read gtf
gtf_df <- read_delim(file = gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# read gtf info
gtf_info <- read_csv(file = gtf_info_path)

# open connection to .fasta file
fa <- open(FaFile(fasta_path))

# read inton sequences
sequences <- Biostrings::readDNAStringSet(filepath = seq_path)

######################################################## MAIN CODE
### get exons from .gtf, take only protein coding genes
gtf <- gtfToGRanges(gtf_df, filter = "exon")
gtf <- gtf[gtf$gene_biotype == "protein_coding"]
# gtf <- gtf[gtf$gene_id == "ENSMUSG00000036698"]

### takes one transcript for each gene (the one with the most exons)
# get unique gene_id-transcript_id combinations
gids <- unique(values(gtf)[c("gene_id", "transcript_id")])

# split exons to transcripts, order by number of exons
gtf_trans <- 
  split(gtf, gtf$transcript_id) %>% 
  .[order(-elementNROWS(.))]

# order unique IDs
gids <- gids[match(names(gtf_trans), gids$transcript_id), ]

# take only those transcripts with most exons
gtf_trans <- gtf_trans[!duplicated(gids$gene_id)]


### enumerate exons in transcripts
## transcripts with single exon
gtf_trans_single <- 
  gtf_trans[elementNROWS(gtf_trans) == 1] %>% 
  unlist(.)
values(gtf_trans_single) <- gtf_trans_single$transcript_id
names(mcols(gtf_trans_single)) <- "transcript_id"
gtf_trans_single$ex.num <- 1
gtf_trans_single$ex.tot <- 1
gtf_trans_single <- split(gtf_trans_single, gtf_trans_single$transcript_id)

## transcripts with more than one exon
gtf_trans <- 
  gtf_trans[elementNROWS(gtf_trans) > 1] %>% 
  reduce(., min.gapwidth = 100) %>% 
  unlist(.)
gtf_trans$transcript_id <- names(gtf_trans)
names(gtf_trans) <- NULL

# create data.frame with enumerated exons
values_df <- 
  values(gtf_trans) %>% 
  as.data.frame(.) %>% 
  as.tibble(.) %>% 
  dplyr::mutate(strand = as.character(strand(gtf_trans))) %>% 
  dplyr::group_by(transcript_id) %>% 
  mutate(count = n(), 
         idx = ifelse(strand == "+", 1:n(), n():1)) %>% 
  dplyr::ungroup(.)

# add enumeration to transcripts
gtf_trans$ex.num <- values_df$idx
gtf_trans$ex.tot <- values_df$count
gtf_trans <- split(gtf_trans, gtf_trans$transcript_id)


### get positions of splice donors and acceptors
# get individual exons
splice <- unlist(gtf_trans)

# # splice donor - 5 nt in exon + 15 nt downstream
# splice_don <- GenomicRanges::resize(splice, fix = "end", width = 5)
# splice_don <- splice_don[splice_don$ex.num != splice_don$ex.tot]
# splice_don <- resize(splice_don, width = 5, fix = "start")
# splice_don <- unique(splice_don)
# splice_don <- splice_don[countOverlaps(splice_don, splice_don) == 1]
# splice_don_seq <- Rsamtools::getSeq(x = fa, splice_don)

# splice acceptor - 5 nt in exon + 10 nt upstream
splice_acc <- GenomicRanges::resize(splice, fix = "start", width = 2)
splice_acc <- splice_acc[splice_acc$ex.num != 1]	
splice_acc <- resize(splice_acc, width = 15, fix = "end")
splice_acc <- unique(splice_acc)
splice_acc <- splice_acc[countOverlaps(splice_acc, splice_acc) == 1]
# splice_acc_seq <- Rsamtools::getSeq(x = fa, splice_acc)
# Biostrings::writeXStringSet(x = splice_acc_seq, filepath = file.path(outpath, "all.splice_acceptors.fa"))
# Biostrings::writeXStringSet(x = splice_acc_seq, filepath = file.path(outpath, "Ago2.splice_acceptors.fa"))
splice_acc_seq <- Biostrings::readDNAStringSet(filepath = file.path(outpath, "Ago2.splice_acceptors.fa"))

### splice acceptor logo
# make logo
sample_pcm <- 
  as.matrix(splice_acc_seq) %>% 
  apply(., 2, str_c, collapse = "") %>% 
  Biostrings::DNAStringSet(.) %>% 
  Biostrings::letterFrequency(., letters = c("A", "C", "G", "T")) %>% 
  apply(., 1, proportion) %>% 
  seqLogo::makePWM(.)

# plot logo
png(filename = file.path(outpath, "mmusculus.splice_acceptor.Ago2.scaled.logo.png"), width = 1500, height = 400)
seqLogo::seqLogo(sample_pcm, ic.scale = T)
dev.off()

### match splice acceptor PWM with intron3_4
# filter N and other ambiguous characters
alph_freq <- Biostrings::alphabetFrequency(splice_acc_seq)
splice_acc_seq <- splice_acc_seq[rowSums(alph_freq[, 5:ncol(alph_freq)]) == 0, ]

# create PWM
splice_acc_pwm <- Biostrings::PWM(x = splice_acc_seq, type = "log2probratio")

# match PWM
splice_acc_ranges <- 
  Biostrings::matchPWM(pwm = splice_acc_pwm, subject = sequences[[12]]) %>%
  # matchPattern(pattern = "AG", subject = sequences[[12]], fixed = T) %>% 
  ranges(.) %>% 
  reduce(.) %>% 
  as.data.frame(.) %>% 
  as.tibble(.) %>% 
  dplyr::mutate(match = "splice_acceptor")


### find branchpoint pattern = YTRAC/CTRAY
branchpoint_ranges <- 
  matchPattern(pattern = "YTRAC", subject = sequences[[12]], fixed = F) %>% 
  ranges(.) %>% 
  reduce(.) %>% 
  as.data.frame(.) %>% 
  as.tibble(.) %>% 
  dplyr::mutate(match = "branchpoint")
  

### bind all, add reference range, get absolute positions, create GRangesList
# create reference coordinates
reference_range <- GenomicRanges::GRanges(seqnames = "chr15", 
                                          ranges = IRanges(start = 73131042, end = 73137993), 
                                          strand = "-", 
                                          match = "reference")

# bind to one df
all_ranges <- 
  rbind(splice_acc_ranges, branchpoint_ranges) %>% 
  dplyr::rename(start_rel = start, end_rel = end) %>% 
  dplyr::mutate(strand = "-", 
                start = start(reference_range), 
                end = end(reference_range)) %>% 
  mutate_cond(strand == "+", start_abs = start + start_rel - 1, end_abs = start_abs + width - 1) %>% 
  mutate_cond(strand == "-", start_abs = end - end_rel + 1 , end_abs = start_abs + width - 1) %>% 
  dplyr::select(start = start_abs, end = end_abs, match, strand) %>% 
  dplyr::mutate(seqnames = "chr15") %>% 
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T) %>% 
  GenomicRanges::split(., .$match)
  
  
### plot
# set different ranges
Ago2_range <- GRanges(seqnames = "chr15", IRanges(start = 73101625, end = 73184947), strand = "-")
intron3_4_range <- GRanges(seqnames = "chr15", IRanges(start = 73131042, end = 73137993), strand = "-")
extra_exon_range <- GRanges(seqnames = "chr15", IRanges(start = 73134498, end = 73134857 + 200))
splice_acceptor_range <- GRanges(seqnames = "chr15", IRanges(start = 73134849, end = 73134863))

# prepare tracks
Fugaku_exp <- autoplot(bam_path, which = reference_range, geom = "line")
Fugaku_exp2 <- autoplot(bam_path, which = reference_range, geom = "gapped.pair")
splice_acc <- autoplot(all_ranges$splice_acceptor, which = reference_range, color = "red")
branchpoint <- autoplot(all_ranges$branchpoint, which = reference_range, color = "blue")
txdb <- autoplot(TxDb.Mmusculus.UCSC.mm10.ensGene, which = reference_range)
genome <- autoplot(BSgenome.Mmusculus.UCSC.mm10, which = reference_range)

# save plot
png(filename = file.path(outpath, "ggbio.splice_acceptor.Ago2.super_small.png"), width = 1000, height = 800)
ggbio::tracks(`Fugaku GV coverage` = Fugaku_exp,
              `Fugaku GV reads` = Fugaku_exp2,
              `Splice acceptor` = splice_acc,
              branchpoint = branchpoint,
              reference = txdb,
              genome = genome,
              heights = c(0.5, 0.3, 0.1, 0.1, 0.1, 0.1), 
              theme = theme_clear(), 
              label.text.cex = 0.75, 
              xlim = splice_acceptor_range)
dev.off()

