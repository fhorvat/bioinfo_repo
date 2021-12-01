### INFO: 
### DATE: Sun Nov 24 12:22:19 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/siRNA_pseudogenes")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(purrr)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(DESeq2)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/Yang_2016_SciAdv_GSE83581/Data/Mapped/STAR_mm10.4_mismatches"

# set outpath
outpath <- getwd()

# bam paths
bam_path <- list.files(inpath, ".*bam$", full.names = T)

# fasta paths
fasta_path <- list.files(file.path(inpath, "transcript_index"), ".*mRNA.fa", full.names = T)

######################################################## READ DATA
# read bam
bam_gr <- purrr::map(bam_path, function(path){
  
  # get read sequences, filter 
  bam_galign <- 
    GenomicAlignments::readGAlignments(file = path, param = ScanBamParam(tag = c("nM", "NH"), what = c("qname", "seq"))) %>% 
    .[str_detect(cigar(.), "2[1-3]M") & !(str_detect(cigar(.), "I|D"))]
  
}) %>%
  do.call(c, .)

# read fasta
genomic_fasta <- Biostrings::readDNAStringSet(filepath = fasta_path)
names(genomic_fasta) <- str_remove(names(genomic_fasta), " .*")

######################################################## MAIN CODE
# get all read names
read_names <- 
  mcols(bam_gr)$qname %>% 
  unique(.)

# get all reads which are sense to Sirena1 (=antisense to Gm31941)
bam_gm <- bam_gr[seqnames(bam_gr) == "Gm31941" & strand(bam_gr) == "-"]

# get all the reads which are sense to either Elob 
bam_elob_elobl <- bam_gr[seqnames(bam_gr) %in% c("Elob", "Elobl") & strand(bam_gr) == "+"]


### get mismatch counts
purrr::map(c("Gm31941", "Elob_Elobl"), function(gene_name){
  
  # which sequence
  if(gene_name == "Gm31941"){
    bam_galign <- bam_gr[seqnames(bam_gr) == gene_name & strand(bam_gr) == "-"]
  }else{
    bam_galign <- bam_gr[seqnames(bam_gr) %in% c("Elob", "Elobl") & strand(bam_gr) == "+"]
  }
  
  # get sequences
  bam_seq <- GenomicAlignments::sequenceLayer(x = mcols(bam_galign)$seq, cigar = cigar(bam_galign))
  
  # get genomic ranges
  genomic_coords <- 
    GRanges(bam_galign) %>% 
    unstrand(.)
  
  # get genomic sequences
  genomic_seq <- genomic_fasta[genomic_coords]
  
  # align sequences
  bam_genomic_align <- pairwiseAlignment(bam_seq, genomic_seq)
  
  # get mismatch table
  alignment_mismatch <- mismatchTable(bam_genomic_align)
  
  # filter and tidy mismatch table
  alignment_mismatch_tb <- 
    alignment_mismatch %>% 
    as_tibble(.) %>% 
    dplyr::mutate(read_width = width(bam_gm)[PatternId]) %>% 
    dplyr::select(read_base = PatternSubstring, genome_base = SubjectSubstring, mismatch_position = PatternStart, read_width) %>% 
    dplyr::filter(!(mismatch_position == read_width & read_base == "A")) %>% 
    dplyr::count(read_base, genome_base) %>% 
    dplyr::filter_at(vars(read_base, genome_base), all_vars(. != "N")) %>% 
    tidyr::pivot_wider(id_cols = read_base, names_from = genome_base, values_from = n) %>%
    dplyr::select(read_base, A, C, G, `T`) %T>% 
    readr::write_csv(., file.path(outpath, str_c(gene_name, ".sense_siRNA.mismatches.csv")))
  
  # sum total number of mismatches in all reads
  alignment_mismatch_tb %>% 
    pivot_longer(-read_base, names_to = "genome_base", values_to = "count", values_drop_na = T) %$%
    count %>% 
    sum(.)
    
})

