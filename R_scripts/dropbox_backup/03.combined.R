#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other) in chunks
### DATE: Mon Mar 04 14:51:59 2019
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd(".")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(data.table)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
## counts reads in each category 
classReads <- function(bam, yield = 1000000, isFirstInPair = NA){
  
  # initialize vector to hold number of reads in each category
  count_sums <- c("rRNA" = 0, "repeats" = 0, "exon" = 0, "other" = 0)
  
  # initialize empty data.tabe
  read_counts_empty <- data.table(hits = 1:3)
  
  # open connection to bam file in chunks
  bamfile <- BamFile(bam, yieldSize = yield)
  open(bamfile)
  
  # load chunks of reads from bam file and classify each alignment
  while(length(chunk <- readGAlignmentsList(bamfile, param = ScanBamParam(what = "qname", flag = scanBamFlag(isFirstMateRead = isFirstInPair, isSecondaryAlignment = F))))) {
    
    # unlist, set names of read
    chunk <- unlist(chunk)
    names(chunk) <- mcols(chunk)$qname
    
    # get number of unique reads in chunk
    unique_reads_number <- length(names(chunk))
    
    # transform to grglist (which gets ranges of only alignment part), unlist to GRanges
    chunk <-
      GenomicRanges::grglist(chunk) %>%
      unlist(.)
    
    # count reads overlaping rDNA, repeats and exons
    read_counts <- findOverlaps(chunk, features_list, ignore.strand = T, type = "within")
    
    # create data.table with counts
    read_counts <-
      data.table(hits = subjectHits(read_counts),
                 reads = names(chunk[queryHits(read_counts)])) %>%
      .[, list(hits = min(hits)), by = "reads"] %>%
      .[order(hits), list(n = .N), by = "hits"] %>%
      .[read_counts_empty, on = "hits"] %>%
      .[is.na(n), n := 0] %$%
      n
    
    # calculate how many reads are not overlaping any class
    read_counts_other <- unique_reads_number - sum(read_counts)
    
    # add other reads to read counts
    read_counts <- c(read_counts, read_counts_other)
    
    # add to count sums
    count_sums <- count_sums + read_counts
    
  }
  
  # close connection to .bam
  close(bamfile)
  
  # return vector with read names
  return(count_sums)
  
}

######################################################## PATH VARIABLES
# set outpath, get paths for gtf and bam files
outpath <- getwd()

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# arguments from command line
single_end <- as.logical(args$single_end)
bam_path <- args$bam_path
bam_name <- args$bam_name
experiment_name <- args$experiment_name
features_exons <- args$features_exons
features_rmsk <- args$features_rmsk
features_geneInfo <- args$features_geneInfo
features_mirbase <- args$features_mirbase
class_algorithm <- args$class_algorithm

single_end='TRUE'
threads=''
bam_path='./s_12oocytes_Mov10l1_WT_So834_12xoo_r1.SE.bam'
bam_name='s_12oocytes_Mov10l1_WT_So834_12xoo_r1.SE'
features_exons='/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/ensembl.99.MesAur1.0.20200415.Siomi.UCSCseqnames.reducedExons.RDS'
features_rmsk='/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/rmsk.Siomi.20200701.clean.fa.out.gz'
features_geneInfo='/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/ensembl.99.MesAur1.0.20200415.Siomi.UCSCseqnames.geneInfo.csv'
features_mirbase=''
class_algorithm='loop'


######################################################## READ DATA
### set feature classes
class_hier <- c("miRNA.mature.sense", "miRNA.other.sense",
                "protein_coding.sense",
                "tRNA", "rRNA", 
                "SINE", "LINE", "LTR", "other_repeat",
                "annotated_pseudogene",
                "other", "not_annotated")

# read gene info
genes_info <- readr::read_csv(features_geneInfo)

# read exons
exons_gr <-
  readRDS(file = features_exons) %>%
  tibble::as_tibble(.) %>%
  dplyr::left_join(., genes_info %>% dplyr::select(gene_id, gene_biotype), by = "gene_id") %>%
  dplyr::select(seqnames:strand, gene_id, gene_biotype) %>%
  dplyr::mutate(gene_biotype =
                  replace(gene_biotype, gene_biotype == "Mt_rRNA", "rRNA") %>%
                  replace(., str_detect(., "pseudogene"), "annotated_pseudogene") %>%
                  replace(., . == "miRNA", "miRNA.other")) %>%
  GenomicRanges::GRanges(.)

# read repeatMasker
rmsk_gr <-
  read_delim(file = features_rmsk, delim = "\t", col_types = cols(start = col_double(), end = col_double())) %>%
  dplyr::mutate(gene_id = str_c(seqnames, ":", start, "-", end, "|", strand, "|", repName)) %>%
  dplyr::select(seqnames:strand, gene_id, gene_biotype = repClass) %>%
  dplyr::mutate(gene_biotype = replace(x = gene_biotype,
                                       list = !(gene_biotype %in% c("rRNA", "tRNA", "LINE", "SINE", "LTR")),
                                       values = "other_repeat")) %>%
  GenomicRanges::GRanges(.)

# read miRBase gff or get miRNA coordinates from ENSEMBL
if(length(features_mirbase) > 0){
  
  # read and clean miRBase miRNA annoation
  mirna_gr <- rtracklayer::import.gff(con = features_mirbase)
  mirna_gr <- mirna_gr[mcols(mirna_gr)$type == "miRNA"]
  mcols(mirna_gr) <- mcols(mirna_gr)[, c("ID")]
  names(mcols(mirna_gr)) <- "gene_id"
  mcols(mirna_gr)$gene_biotype <- "miRNA.mature"
  
}else{
  
  # clean ENSEMBL miRNA annoation
  mirna_gr <- exons_gr[mcols(exons_gr)$gene_biotype == "miRNA.other"]
  exons_gr <- exons_gr[mcols(exons_gr)$gene_biotype != "miRNA.other"]
  
}

# join repeatMasker, exons and miRNA annotation
features_gr <- c(rmsk_gr, exons_gr, mirna_gr)

# set all other categorie to "other"
mcols(features_gr)$gene_biotype <- replace(mcols(features_gr)$gene_biotype, 
                                           !(mcols(features_gr)$gene_biotype) %in% class_hier, 
                                           "other")

# join to one list
features_list <- split(features_gr, mcols(features_gr)$gene_biotype)

######################################################## MAIN CODE


# classify reads in bam file
reads_class_counts <- classReads(bam = bam_path, yield = 1000000, isFirstInPair = ifelse(!single_end, T, NA))

# read second reads in pair if experiment was paired end
if(!single_end){
  
  # class reads
  reads_class_counts_second <- classReads(bam = bam_path, yield = 1000000, isFirstInPair = F)
  
  # sum with counts
  reads_class_counts <- reads_class_counts + reads_class_counts_second
  
  # get number of FRAGMENTS which are mapping to genome except rDNA (for normalizing)
  genome.mapped_minus_rDNA <-
    (sum(reads_class_counts) - unname(reads_class_counts["rRNA"])) %>%
    magrittr::divide_by(., 2) %>%
    round(.)
  
  # add number of reads in each class to table, save
  reads_class_sum_final <-
    tibble(read_group = names(reads_class_counts), count = reads_class_counts) %>%
    tidyr::spread(., read_group, count) %>%
    dplyr::mutate(sample_id = bam_name,
                  total = rRNA + repeats + exon + other,
                  genome.mapped_minus_rDNA = genome.mapped_minus_rDNA) %>%
    dplyr::select(sample_id, rRNA, repeats, exon, other, total, genome.mapped_minus_rDNA) %T>%
    readr::write_delim(., file = file.path(outpath, str_c(bam_name, ".read_stats.txt")), delim = "\t")
  
}else{
  
  # add number of reads in each class to table, save
  reads_class_sum_final <-
    tibble(read_group = names(reads_class_counts), count = reads_class_counts) %>%
    tidyr::spread(., read_group, count) %>%
    dplyr::mutate(sample_id = bam_name,
                  total = rRNA + repeats + exon + other,
                  genome.mapped_minus_rDNA = total - rRNA) %>%
    dplyr::select(sample_id, rRNA, repeats, exon, other, total, genome.mapped_minus_rDNA) %T>%
    readr::write_delim(., file = file.path(outpath, str_c(bam_name, ".read_stats.txt")), delim = "\t")
  
}
