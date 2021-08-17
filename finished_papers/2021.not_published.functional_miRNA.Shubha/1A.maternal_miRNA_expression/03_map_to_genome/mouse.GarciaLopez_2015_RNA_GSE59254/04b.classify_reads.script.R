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
classReadsFactor <- function(bam_path, subject_ranges, yield = 1000000, isFirstInPair = NA){
  
  # get number of alignments in bam file
  read_number <-
    Rsamtools::countBam(file = bam_path, param = ScanBamParam(flag = scanBamFlag(isFirstMateRead = isFirstInPair))) %$%
    records
  
  # initialize alignment class vector
  align_class <- rep("placeholder", read_number)
  names(align_class) <- rep("nameplaceholder", read_number)
  
  # open connection to bam file in chunks
  bamfile <- BamFile(bam_path, yieldSize = yield)
  open(bamfile)
  
  # load chunks of alignments from bam file and classify each alignment
  while(length(chunk <- readGAlignmentsList(bamfile, param = ScanBamParam(what = "qname", flag = scanBamFlag(isFirstMateRead = isFirstInPair))))) {
    
    # unlist, set names of reads, transform to grglist (which gets ranges of only alignment part)
    chunk <- unlist(chunk)
    names(chunk) <- mcols(chunk)$qname
    chunk <-
      GenomicRanges::grglist(chunk) %>%
      unlist(.)
    
    # find overlaps between two GRanges - chunk of alignments and annotation
    hits <- findOverlaps(chunk, subject_ranges, ignore.strand = T)
    
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

### class reads and save table
reads_class <- classReadsFactor(bam_path = bam_path,
                                subject_ranges = features_gr,
                                yield = 1000000,
                                isFirstInPair = NA)

# # create class hierarchy table
# hier_dt <- data.table(position = 1:length(class_hier), 
#                       read_group = class_hier)
# data.table::setkey(hier_dt, "read_group")
# 
# # create class table
# read_class_dt <- data.table(read_name = reads_class, 
#                             read_group = names(reads_class))
# data.table::setkey(read_class_dt, "read_group")
# 
# # get highest hit for each read
# class_sum <- 
#   read_class_dt[hier_dt] %>% 
#   .[, list(position = min(position)), by = "read_name"] %>% 
#   .[order(position), list(count = .N), by = "position"] %>%
#   .[hier_dt, on = "position"] %>%
#   .[is.na(count), count := 0] %>% 
#   .[, `:=`(sample_id = bam_name)] %>%
#   .[order(position), ] %>% 
#   .[, position := NULL] %>% 
#   setcolorder(., c("read_group", "count", "sample_id"))
# 
# # write
# readr::write_delim(x = class_sum, path = file.path(outpath, str_c(bam_name, "read_stats", "txt", sep = ".")), delim = "\t")


### loop through unique names and sum read categories
# get unique names
reads_class_unique <- unique(reads_class)

# create sum table
reads_class_sum <- tibble(read_group = class_hier,
                          count = 0)
 
# set loop step
loop_step <- 100000

# set loop index
loop_index <- seq(0, length(reads_class_unique), by = loop_step)
if(loop_index[length(loop_index)] < length(reads_class_unique)){
  loop_index <- c(loop_index, length(reads_class_unique))
}

# loop
for(n in 1:(length(loop_index) - 1)){
  
  # subset vector
  reads_class_subset <- reads_class[reads_class %in% reads_class_unique[(loop_index[n] + 1):loop_index[n + 1]]]
  
  # summarize read classification
  read_class_vector <-
    factor(x = names(reads_class_subset), levels = class_hier) %>%
    as.numeric(.) %>%
    split(., factor(reads_class_subset)) %>%
    purrr::map(., min) %>%
    unlist(.)
  
  # summarize
  reads_class_sum_subset <-
    tibble::tibble(read_id = names(read_class_vector),
                   read_group = class_hier[read_class_vector]) %>%
    unique(.) %>%
    dplyr::group_by(read_group) %>%
    dplyr::summarise(count = n())
  
  # join with sum vector
  reads_class_sum %<>%
    dplyr::left_join(., reads_class_sum_subset, by = "read_group") %>%
    dplyr::mutate(count.y = replace(count.y, is.na(count.y), 0),
                  count = count.x + count.y) %>%
    dplyr::select(read_group, count)
  
}
 
# add info about sample and experiment, save
reads_class_sum %<>%
  dplyr::mutate(sample_id = bam_name,
                experiment = experiment_name) %T>%
  readr::write_delim(x = ., path = file.path(outpath, str_c(bam_name, "read_stats", "txt", sep = ".")), delim = "\t")
