#!/home/students/fhorvat/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other)
### DATE: 03. 09. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("%OUTPATH")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tibble)
library(data.table)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)

######################################################## PATH VARIABLES
gtf_path <- "/common/WORK/fhorvat/reference/golden_hamster/mesAur1/GCF_000349665.1_MesAur1.0_genomic.gtf"
repeatmasker_path <- "/common/WORK/fhorvat/reference/golden_hamster/mesAur1/mesAur1.RepeatMasker.Refseq.txt"
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"

outpath <- getwd()
file <- "%FILE"
file_name <- str_replace_all(basename(file), ".Aligned.sortedByCoord.out.bam", "")

######################################################## SOURCE FILES
source(file.path(lib_path, "GffToGRanges.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
# refSeq gtf
gtf_exons <-
  read_delim(file = gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c")) %>%
  GffToGRanges(., "exon")

# repeatMasker
rmsk_gr <-
  read_delim(file = repeatmasker_path, delim = "\t", col_names = F) %>%
  magrittr::set_colnames(c("seqnames", "start", "end", "strand", "repName", "repClass_repFamily")) %>% 
  tidyr::separate(repClass_repFamily, into = c("repClass", "repFamily"), sep = "/") %>% 
  dplyr::mutate(strand = replace(strand, strand == "C", "*")) %>% 
  # dplyr::mutate(full_pos = str_c(seqnames, ":", start, "-", end, "|", strand, "|", repClass, "|", repName)) %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

# repeatMasker - rRNA
rmsk_rRNA <- rmsk_gr[rmsk_gr$repClass == "rRNA"]

######################################################## MAIN CODE
### count reads over features hierarchically
# complete match: rDNA -> repeat -> exon -> other
# read in the .bam, set unique names for each mate in fragment
# readGAlignmentPairs(file, param = ScanBamParam(which = GRanges("NW_004801604.1", IRanges(1, 1500000))), use.names = TRUE) 
  
gbam <-
  readGAlignmentPairs(file, use.names = TRUE) %>%
  grglist(.) %>%
  unlist(.)
names(gbam) <- str_c(names(gbam), 1:2)

# count reads in each category (rDNA -> repeat -> exon -> other) with subsequent filtering
ordered_filters <-
  list(rmsk_rRNA, rmsk_gr, gtf_exons) %>%
  magrittr::set_names(c("rRNA", "repeats", "exon"))
reads_sum <- tibble(sample = file_name)

for(overlap_subject in names(ordered_filters)){
  
  # overlap reads
  overlap_reads <-
    IRanges::subsetByOverlaps(query = gbam, subject = ordered_filters[[overlap_subject]], type = "within", ignore.strand = T) %>%
    names(.) %>%
    unique(.)
  
  # add columnt to sum table
  reads_sum %<>%
    mutate(overlap_subject = length(overlap_reads)) %>%
    data.table::setnames(., old = ncol(.), new = overlap_subject)
  
  # filter reads out
  gbam <- gbam[!(names(gbam) %in% overlap_reads)]
  
}

# set everything else to "other" and sum total
reads_sum %<>%
  dplyr::mutate(other = length(unique(names(gbam)))) %>%
  dplyr::mutate(total = rowSums(.[, -1])) %T>%
  readr::write_delim(., path = file.path(outpath, str_c(file_name, "_read_stats.txt")), delim = "\t")
