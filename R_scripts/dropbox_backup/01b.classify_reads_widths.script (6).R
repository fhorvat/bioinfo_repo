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
classReadsWidth <- function(bam_path, subject_ranges, yield = 1000000, isFirstInPair = NA){
  
  # get number of alignments in bam file
  read_number <-
    Rsamtools::countBam(file = bam_path, param = ScanBamParam(flag = scanBamFlag(isFirstMateRead = isFirstInPair,
                                                                                 isSecondaryAlignment = F))) %$%
    records
  
  # create class hierarchy table
  hier_dt <- data.table(position = 1:length(class_hier),
                        read_group = class_hier)
  data.table::setkey(hier_dt, "read_group")
  
  # get table with all combination of positions and widths
  read_width_sum <-
    expand.grid(1:(length(class_hier) + 1), 19:32) %>%
    as.data.table(.) %>%
    setnames(., c("position", "read_width")) %>%
    .[, `:=`(count = 0)] %>%
    .[]
  data.table::setkey(read_width_sum, position, read_width)
  
  # open connection to bam file in chunks
  bamfile <- BamFile(bam_path, yieldSize = yield)
  open(bamfile)
  
  # load chunks of alignments from bam file and classify each alignment
  while(length(chunk <- readGAlignmentsList(bamfile,
                                            param = ScanBamParam(what = "qname",
                                                                 flag = scanBamFlag(isFirstMateRead = isFirstInPair,
                                                                                    isSecondaryAlignment = F))))) {
    
    # unlist, set names of reads, transform to grglist (which gets ranges of only alignment part)
    chunk <- unlist(chunk)
    names(chunk) <- mcols(chunk)$qname
    chunk <-
      GenomicRanges::grglist(chunk) %>%
      unlist(.)
    
    # find overlaps between two GRanges - chunk of alignments and annotation. Get only sense alignments
    hits <- findOverlaps(chunk, subject_ranges, ignore.strand = F)
    
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
      
      # create class table
      read_class_dt <- data.table(read_name = names(read_hits),
                                  read_group = mcols(subject_hits)$gene_biotype_precise,
                                  read_width = width(read_hits))
      data.table::setkey(read_class_dt, "read_group")
      
      # get highest hit for each read, sum per read width and hierarchy
      class_sum <-
        read_class_dt[hier_dt] %>%
        .[!is.na(read_name), ] %>%
        .[, list(position = min(position), read_width = min(read_width)), by = "read_name"] %>%
        .[, list(count = .N), by = c("position", "read_width")] %>%
        .[hier_dt, on = "position"] %>%
        .[is.na(count), count := 0]
      
    }else{
      
      # create empty data.table
      class_sum <- data.table()
      
    }
    
    
    # get reads which were mapped, but not annotated
    if(!(all(names(chunk) %in% names(read_hits)))){
      
      # create data.table with widths of reads
      reads_not_annotated <-
        data.table(read_width = width(chunk)[!(names(chunk) %in% names(read_hits))]) %>%
        .[, list(count = .N), by = c("read_width")] %>%
        .[, `:=`(read_group = "not annotated",
                 position = (nrow(hier_dt) + 1))] %>%
        setcolorder(., c("position", "read_width", "count", "read_group")) %>%
        .[]
      
    }else{
      
      # create empty data.table
      reads_not_annotated <- data.table()
      
    }
    
    # bind
    class_sum <- rbind(class_sum, reads_not_annotated)
    data.table::setkey(class_sum, position, read_width)
    
    # join with table from previous chunk
    read_width_sum[class_sum, count := count + i.count]
    
  }
  
  # close connection to .bam
  close(bamfile)
  
  # add names
  read_width_sum[hier_dt, read_group := i.read_group, on = "position"]
  read_width_sum[is.na(read_group), `:=`(read_group = "not annotated")]
  
  # check if we counted all the reads
  if(read_number != sum(read_width_sum$count)) stop("You missed some reads mate!")
  
  # return vector with read names
  return(as_tibble(read_width_sum))
  
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
features_ensembl <- args$features_ensembl
features_rmsk <- args$features_rmsk

######################################################## READ DATA
# read .gtf
gtf_gr <- rtracklayer::import(features_ensembl)

# read repeatMasker
rmsk_tb <- read_delim(file = features_rmsk, delim = "\t", col_types = cols(start = col_double(), end = col_double()))

######################################################## MAIN CODE
### prepare annotation
# repeatMasker table
rmsk_key_tb <-
  rmsk_tb %>%
  dplyr::select(repName, repClass, repFamily) %>%
  dplyr::distinct(.) %>%
  dplyr::mutate(gene_biotype = NA,
                gene_biotype = replace(gene_biotype, repClass == "rRNA", "rRNA"),
                gene_biotype = replace(gene_biotype, repClass == "tRNA", "tRNA"),
                gene_biotype = replace(gene_biotype, repClass == "LINE", "LINE"),
                gene_biotype = replace(gene_biotype, repClass == "SINE", "SINE"),
                gene_biotype = replace(gene_biotype, repClass == "LTR", "LTR"),
                gene_biotype = replace(gene_biotype, str_detect(repName, "^IAP"), "IAP"),
                gene_biotype = replace(gene_biotype, repClass %in% c("DNA", "DNA?", "Simple_repeat"), "DNA/simple repeat"),
                gene_biotype = replace(gene_biotype, repClass %in% c("snRNA", "srpRNA", "scRNA", "RNA"), "non-coding RNA RepeatMasker"),
                gene_biotype = replace(gene_biotype, is.na(gene_biotype), "other annotated RepeatMasker")) %>%
  dplyr::mutate(gene_biotype_precise = NA,
                gene_biotype_precise = ifelse(gene_biotype %in% c("LINE", "SINE", "LTR"), repFamily, NA),
                gene_biotype_precise = replace(gene_biotype_precise, (gene_biotype == "SINE" & gene_biotype_precise == "tRNA"), "SINE tRNA"),
                gene_biotype_precise = replace(gene_biotype_precise,
                                               !(gene_biotype_precise %in% c("L1", "L2", "CR1", "RTE-X",
                                                                             "ERVK", "ERVL-MaLR", "ERVL", "ERV1", "LTR", "Gypsy",
                                                                             "Alu", "B2", "B4", "MIR", "ID", "SINE tRNA", "Deu")),
                                               NA)) %>%
  dplyr::mutate(gene_biotype_precise = replace(gene_biotype_precise,
                                               ((gene_biotype == "LINE") & (is.na(gene_biotype_precise))),
                                               "LINE other"))

# clean repeatMasker
rmsk_gr <-
  rmsk_tb %>%
  dplyr::mutate(gene_id = str_c(seqnames, ":", start, "-", end, "|", strand, "|", repName)) %>%
  dplyr::left_join(., rmsk_key_tb, by = "repName") %>%
  dplyr::select(seqnames:strand, gene_id, gene_biotype, gene_biotype_precise) %>%
  GenomicRanges::GRanges(.)

# get miRNA annotation
mirna_gr <- gtf_gr[mcols(gtf_gr)$gene_biotype == "miRNA"]
mirna_gr <- mirna_gr[mcols(mirna_gr)$type == "gene"]
mcols(mirna_gr)$gene_biotype_precise <- NA

# clean ENSEMBL
ensembl_gr <-
  gtf_gr %>%
  as_tibble(.) %>%
  dplyr::filter(gene_biotype != "miRNA", type == "exon") %>%
  dplyr::mutate(gene_biotype = replace(gene_biotype, gene_biotype == "Mt_rRNA", "rRNA")) %>%
  dplyr::mutate(gene_biotype = replace(gene_biotype, gene_biotype == "Mt_tRNA", "tRNA")) %>%
  dplyr::mutate(gene_biotype = replace(gene_biotype,
                                       gene_biotype %in% c("snRNA", "lncRNA", "snoRNA", "miscRNA",
                                                           "ribozyme", "scaRNA", "sRNA", "scRNA",
                                                           "Mt_tRNA"),
                                       "non-coding RNA ENSEMBL"),
                gene_biotype = replace(gene_biotype,
                                       !(gene_biotype %in% c("rRNA", "non-coding RNA ENSEMBL", "protein_coding")),
                                       "other annotated ENSEMBL"),
                gene_biotype = replace(gene_biotype, gene_biotype == "protein_coding", "protein coding ENSEMBL"),
                gene_biotype_precise = NA) %>%
  GRanges(.)

# join miRNA, repeatMasker and ENSEMBL annotation
features_gr <- c(mirna_gr, rmsk_gr, ensembl_gr)
mcols(features_gr) <- mcols(features_gr)[, c("gene_id", "gene_biotype", "gene_biotype_precise")]
mcols(features_gr)$gene_biotype_precise <- ifelse(is.na(mcols(features_gr)$gene_biotype_precise),
                                                  mcols(features_gr)$gene_biotype,
                                                  mcols(features_gr)$gene_biotype_precise)


### set feature classes
class_hier <- c("miRNA",
                "rRNA", "tRNA",
                "L1", "L2", "CR1", "RTE-X", "LINE other",
                "ERVK", "ERVL-MaLR", "ERVL", "ERV1", "LTR", "Gypsy",
                "IAP",
                "Alu", "B2", "B4", "MIR", "ID", "SINE tRNA", "Deu",
                "DNA/simple repeat",
                "non-coding RNA RepeatMasker", "other annotated RepeatMasker",
                "protein coding ENSEMBL", "non-coding RNA ENSEMBL", "other annotated ENSEMBL")

# get features table
features_tb <-
  mcols(features_gr) %>%
  as_tibble(.) %>%
  dplyr::select(read_group = gene_biotype_precise, gene_biotype) %>%
  dplyr::distinct(.)

# class reads
reads_class <- classReadsWidth(bam_path = bam_path,
                               subject_ranges = features_gr,
                               yield = 1000000,
                               isFirstInPair = NA)

# add sample ID to the table, save
reads_class_tb <-
  reads_class %>%
  dplyr::mutate(sample_id = bam_name,
                experiment_name = experiment_name) %>%
  dplyr::left_join(., features_tb, by = "read_group") %T>%
  readr::write_csv(., file = file.path(outpath, str_c(bam_name, "read_class.widths", "csv", sep = ".")))
