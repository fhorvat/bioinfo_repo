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
classReadsHierWidth <- function(bam_path, subject_ranges, yield = 1000000, isFirstInPair = NA, overlap_type = "any"){
  
  # get table with all combination of positions and widths
  read_counts_sum <-
    expand.grid(1:(length(class_hier) + 1), 1:100) %>%
    as.data.table(.) %>%
    setnames(., c("position", "read_width")) %>%
    .[, `:=`(count = 0)] %>%
    .[]
  
  # set key
  data.table::setkey(read_counts_sum, position, read_width)
  
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
        data.table(position = match(gene_biotype_vector, class_hier),
                   reads = names(read_hits), 
                   read_width = width(read_hits)) %>%
        .[, list(position = min(position), read_width = min(read_width)), by = "reads"] %>%
        .[order(position), list(count = .N), by = c("position", "read_width")]
      
    }else{
      
      # create empty data.table
      read_counts <- data.table()
      
    }
    
    # get reads which were mapped, but not annotated
    if(!(all(names(chunk) %in% names(read_hits)))){
      
      # create data.table with widths of reads
      read_counts_not_annotated <-
        data.table(read_width = width(chunk)[!(names(chunk) %in% names(read_hits))]) %>%
        .[, list(count = .N), by = c("read_width")] %>%
        .[, `:=`(position = (length(class_hier) + 1))] %>%
        setcolorder(., c("position", "read_width", "count")) %>%
        .[]
      
    }else{
      
      # create empty data.table
      read_counts_not_annotated <- data.table()
      
    }
    
    # bind
    read_counts <- rbind(read_counts, read_counts_not_annotated)
    data.table::setkey(read_counts, position, read_width)
    
    # join with table from previous chunk
    read_counts_sum[read_counts, count := count + i.count]
    
  }
  
  # close connection to .bam
  close(bamfile)
  
  # return
  return(read_counts_sum)
  
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
features_exons <- args$features_exons
features_rmsk <- args$features_rmsk
features_geneInfo <- args$features_geneInfo
features_mirbase <- args$features_mirbase

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
                "other")


### class reads
reads_class_final <- classReadsHierWidth(bam_path = bam_path,
                                         subject_ranges = features_gr,
                                         yield = 1000000,
                                         isFirstInPair = NA,
                                         overlap_type = "any")

# get number of alignments in bam file
read_number <-
  Rsamtools::countBam(file = bam_path, param = ScanBamParam(flag = scanBamFlag(isFirstMateRead = NA,
                                                                               isSecondaryAlignment = F))) %$%
  records

# check if the counts match the alignment number
if(!(sum(reads_class_final$count) == read_number)){
  
  # stop the script
  stop("You lost some reads mate! Check your script!")
  
}


### create table, save
# create class hierarchy table
class_hier_tb <- tibble(position = 1:(length(class_hier) + 1),
                        read_group = c(class_hier, "not_annotated"))

# add number of reads in each class to table
reads_class_final_sum <- 
  reads_class_final %>% 
  as_tibble(.) %>% 
  dplyr::left_join(., class_hier_tb, by = "position") %>% 
  dplyr::mutate(sample_id = bam_name) %>% 
  dplyr::select(position, read_width, count, read_group, sample_id)

# save table
readr::write_csv(reads_class_final_sum, file = file.path(outpath, str_c(bam_name, "read_class.widths", "csv", sep = ".")))
