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
# class alignments hierarchically
classReadsFactor <- function(bam_path, subject_ranges, yield = 1000000, isFirstInPair = NA, overlap_type = "any"){
  
  # initialize vector to hold number of reads in each category
  count_sums <- rep(0, length(class_hier))
  names(count_sums) <- class_hier
  
  # initialize empty data.tabe
  read_counts_empty <- data.table(hits = 1:length(class_hier))
  
  # open connection to bam file in chunks
  bamfile <- BamFile(bam_path, yieldSize = yield)
  open(bamfile)
  
  # load chunks of alignments from bam file and classify each alignment
  while(length(chunk <- readGAlignmentsList(file = bamfile, 
                                            param = ScanBamParam(what = "qname", 
                                                                 flag = scanBamFlag(isFirstMateRead = isFirstInPair, 
                                                                                    isSecondaryAlignment = F))))) {
    
    # unlist, set names of reads
    chunk <- unlist(chunk)
    names(chunk) <- mcols(chunk)$qname
    
    # get number of unique reads in chunk
    unique_reads_number <- length(names(chunk))
    
    # transform to grglist (which gets ranges of only alignment part)
    chunk <-
      GenomicRanges::grglist(chunk) %>%
      unlist(.)
    
    # find overlaps between two GRanges - chunk of alignments and annotation
    hits <- findOverlaps(chunk, subject_ranges, ignore.strand = T, type = overlap_type)
    
    # class if there are any hits
    if(length(hits) > 0){
      
      # get hits in alignments
      read_hits <-
        extractList(chunk, as(queryHits(hits), "List")) %>%
        unlist(.)
      
      # get hits in annotation
      subject_hits <-
        extractList(subject_ranges, as(subjectHits(hits), "List")) %>%
        unlist(.)
      
      # create vector of gene biotypes
      sense_strand <- ifelse(strand(read_hits) == strand(subject_hits), "sense", "antisense")
      gene_biotype_vector <- subject_hits$gene_biotype
      gene_biotype_vector[str_detect(gene_biotype_vector, "miRNA")] <- str_c(gene_biotype_vector[str_detect(gene_biotype_vector, "miRNA")], ".",
                                                                             sense_strand[str_detect(gene_biotype_vector, "miRNA")])
      gene_biotype_vector[gene_biotype_vector == "protein_coding"] <- str_c(gene_biotype_vector[gene_biotype_vector == "protein_coding"], ".",
                                                                            sense_strand[gene_biotype_vector == "protein_coding"])
      gene_biotype_vector <- replace(gene_biotype_vector,
                                     !(gene_biotype_vector %in% class_hier[class_hier != "other"]),
                                     "other")
      
      # create data.table with counts
      read_counts <-
        data.table(hits = match(gene_biotype_vector, class_hier),
                   reads = names(read_hits)) %>%
        .[, list(hits = min(hits)), by = "reads"] %>%
        .[order(hits), list(n = .N), by = "hits"] %>%
        .[read_counts_empty, on = "hits"] %>%
        .[is.na(n), n := 0] %$%
        n
      
      # calculate how many reads are not overlaping any class
      read_counts_other <- unique_reads_number - sum(read_counts)
      
      # add other reads to read counts
      read_counts[length(read_counts)] <- read_counts_other
      
      # add to count sums
      count_sums <- count_sums + read_counts
      
    }else{
      
      # add all reads in chunk to "not_annotated" category
      read_counts <- rep(0, length(class_hier))
      read_counts[class_hier] <- unique_reads_number
      
      # add to count sums
      count_sums <- count_sums + read_counts
      
    }
    
  }
  
  # close connection to .bam
  close(bamfile)
  
  # return
  return(count_sums)
  
}

# class alignments hierarchically
classReadsFactorOld <- function(bam_path, subject_ranges, yield = 1000000, isFirstInPair = NA, overlap_type = "any"){
  
  # get number of alignments in bam file
  read_number <-
    Rsamtools::countBam(file = bam_path, 
                        param = ScanBamParam(flag = scanBamFlag(isFirstMateRead = isFirstInPair))) %$%
    records
  
  # initialize alignment class vector
  align_class <- rep("placeholder", read_number)
  names(align_class) <- rep("nameplaceholder", read_number)
  
  # open connection to bam file in chunks
  bamfile <- BamFile(bam_path, yieldSize = yield)
  open(bamfile)
  
  # load chunks of alignments from bam file and classify each alignment
  while(length(chunk <- readGAlignmentsList(file = bamfile, 
                                            param = ScanBamParam(what = "qname", 
                                                                 flag = scanBamFlag(isFirstMateRead = isFirstInPair))))) {
    
    # unlist, set names of reads, transform to grglist (which gets ranges of only alignment part)
    chunk <- unlist(chunk)
    names(chunk) <- mcols(chunk)$qname
    chunk <-
      GenomicRanges::grglist(chunk) %>%
      unlist(.)
    
    # find overlaps between two GRanges - chunk of alignments and annotation
    hits <- findOverlaps(chunk, subject_ranges, ignore.strand = T, type = overlap_type)
    
    # class if there are any hits
    if(length(hits) > 0){
      
      # get hits in alignments
      read_hits <-
        extractList(chunk, as(queryHits(hits), "List")) %>%
        unlist(.)
      
      # get hits in annotation
      subject_hits <-
        extractList(subject_ranges, as(subjectHits(hits), "List")) %>%
        unlist(.)
      
      # create vector of gene biotypes
      sense_strand <- ifelse(strand(read_hits) == strand(subject_hits), "sense", "antisense")
      gene_biotype_vector <- subject_hits$gene_biotype
      gene_biotype_vector[str_detect(gene_biotype_vector, "miRNA")] <- str_c(gene_biotype_vector[str_detect(gene_biotype_vector, "miRNA")], ".",
                                                                             sense_strand[str_detect(gene_biotype_vector, "miRNA")])
      gene_biotype_vector[gene_biotype_vector == "protein_coding"] <- str_c(gene_biotype_vector[gene_biotype_vector == "protein_coding"], ".",
                                                                            sense_strand[gene_biotype_vector == "protein_coding"])
      gene_biotype_vector <- replace(gene_biotype_vector,
                                     !(gene_biotype_vector %in% class_hier[class_hier != "other"]),
                                     "other")
      
      # class each unique read ID to category
      read_class_vector <-
        factor(x = gene_biotype_vector,
               levels = class_hier) %>%
        as.numeric(.) %>%
        split(., factor(names(read_hits))) %>%
        purrr::map(., min) %>%
        unlist(.)
      
      # find first element which is not filled
      last_element <-
        which(align_class == "placeholder") %>%
        min(.)
      
      # set range in intialized vector
      class_range <- c(last_element, (last_element - 1 + length(read_class_vector)))
      
      # create named vector (read names and class as name)
      align_class[class_range[1]:class_range[2]] <- names(read_class_vector)
      names(align_class)[class_range[1]:class_range[2]] <- class_hier[read_class_vector]
      
      # filter reads out
      chunk <- chunk[!(names(chunk) %in% names(read_hits))]
      
    }
    
    # set alignements which are not overlaping any category as "not_annotated"
    if(length(chunk) > 0){
      
      # find first element which is not filled
      last_element <-
        which(align_class == "placeholder") %>%
        min(.)
      
      # set range in intialized vector
      class_range <- c(last_element, (last_element - 1 + length(unique(names(chunk)))))
      
      # add to named vector (read names and class as name)
      align_class[class_range[1]:class_range[2]] <- unique(names(chunk))
      names(align_class)[class_range[1]:class_range[2]] <- "not_annotated"
      
    }
    
  }
  
  # close connection to .bam
  close(bamfile)
  
  # remove placeholders from vector
  align_class <- align_class[align_class != "placeholder"]
  
  # return vector with read names
  return(align_class)
  
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
threads='1'
bam_path='./s_BL6Crl_WT_r1.SE.21to23nt.bam'
bam_name='s_BL6Crl_WT_r1.SE.21to23nt'
features_exons='/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.99.GRCm38.p6.20200415.UCSCseqnames.reducedExons.RDS'
features_rmsk='/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/rmsk.mm10.20180919.clean.fa.out.gz'
features_geneInfo='/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.99.GRCm38.p6.20200415.UCSCseqnames.geneInfo.csv'
features_mirbase='/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/miRBase.22.mm10.20181605.gff3'
class_algorithm='loop'

read_number <-
  Rsamtools::countBam(file = bam_path, 
                      param = ScanBamParam(flag = scanBamFlag(isFirstMateRead = NA, 
                                                              isSecondaryAlignment = F))) %$%
  records

######################################################## READ DATA
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
  dplyr::filter(!(repClass %in% c("Simple_repeat", "Low_complexity"))) %>% 
  dplyr::mutate(gene_id = str_c(seqnames, ":", start, "-", end, "|", strand, "|", repName)) %>%
  dplyr::select(seqnames:strand, gene_id, gene_biotype = repClass) %>%
  dplyr::mutate(gene_biotype = replace(x = gene_biotype,
                                       list = !(gene_biotype %in% c("rRNA", "LINE", "SINE", "LTR")),
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

######################################################## MAIN CODE
### set feature classes
class_hier <- c("miRNA.mature.sense", "miRNA.other.sense",
                "protein_coding.sense",
                "rRNA", "SINE", "LINE", "LTR", "other_repeat",
                "annotated_pseudogene",
                "other", "not_annotated")

### class reads
reads_class_final <- classReadsFactor(bam_path = bam_path,
                                      subject_ranges = features_gr,
                                      yield = 1000000,
                                      isFirstInPair = NA,
                                      overlap_type = "any")

### create table, save
# add number of reads in each class to table
reads_class_final_sum <-
  tibble(read_group = names(reads_class_final), 
         count = reads_class_final, 
         sample_id = bam_name)

# # save table
# readr::write_delim(reads_class_final_sum, file = file.path(outpath, str_c(bam_name, ".read_stats.txt")), delim = "\t")


### old algorithm
reads_class <- classReadsFactorOld(bam_path = bam_path,
                                   subject_ranges = features_gr,
                                   yield = 1000000,
                                   isFirstInPair = NA, 
                                   overlap_type = "any")

# create class hierarchy table
hier_dt <- data.table(position = 1:length(class_hier),
                      read_group = class_hier)
data.table::setkey(hier_dt, "read_group")

# create class table
read_class_dt <- data.table(read_name = reads_class,
                            read_group = names(reads_class))
data.table::setkey(read_class_dt, "read_group")

# get highest hit for each read
class_sum <-
  read_class_dt[hier_dt] %>%
  .[, list(position = min(position)), by = "read_name"] %>%
  .[order(position), list(count = .N), by = "position"] %>%
  .[hier_dt, on = "position"] %>%
  .[is.na(count), count := 0] %>%
  .[, `:=`(sample_id = bam_name)] %>%
  .[order(position), ] %>%
  .[, position := NULL] %>%
  setcolorder(., c("read_group", "count", "sample_id"))
