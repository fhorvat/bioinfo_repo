### INFO: annotates repeatMasker and exon overlap of unspliced reads from 1C stage SMRT sequencing
### DATE: 30. 11. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Analysis/unspliced_reads.20171130")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(ggplot2)
library(tibble)
library(data.table)

library(GenomicRanges)

######################################################## PATH VARIABLES
outpath <- getwd()

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "bamToGRangesList.R"))
source(file.path(lib_path, "vfranke", "GffToGRanges.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
### ENSEMBL gtf
# filter, convert to GRanges
gtf_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.gtf.gz"
gtf <- read.table(gtf_path, header = FALSE, sep = "\t", stringsAsFactors = F)
gtf[, 1] <- paste0("chr", gtf[, 1])
gtf <- gtf[!str_detect(gtf[, 1], "NT"),]
gtf[gtf[, 1] == "chrMT", 1] <- "chrM"
gtf <- GffToGRanges(gtf, "exon")

### repeatMasker
rmsk_path <- "/common/WORK/fhorvat/reference/mouse/mm10/UCSC/rmskVIZ.mm10.20171003.clean.OutBaseline.txt.gz" 
rmsk <- 
  readr::read_delim(file = rmsk_path, delim = "\t", col_names = T, col_types = cols(.default = "c")) %>% 
  dplyr::filter(str_detect(repClass, "LTR|LINE|SINE"))
rmsk_gr <- makeGRangesFromDataFrame(rmsk, keep.extra.columns = T)

###  bam files
# bam to GRangesList, make unique names for multimappers
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/SMRT_oocytes_2017/Analysis/unspliced_reads.20171130/unspliced_bams"
bam_list <- list.files(path = bam_path, pattern = "*._1C.*bam$", full.names = T)
reads_grlist <- 
  lapply(bam_list, bamToGRangesList) %>% 
  do.call(what = c, .)
names(reads_grlist) <- make.unique(names(reads_grlist), sep = "_")

######################################################## MAIN CODE
#### get reads info 
# only works on unspliced reads
reads_gr <- unlist(reads_grlist)
reads_gr$read_id <- names(reads_gr)
names(reads_gr) <- NULL
reads_gr <- reads_gr[seqnames(reads_gr) != "chrM"]

# get coordinates of each read in data.frame
reads_df <- 
  reads_gr %>% 
  as.data.frame(.) %>% 
  as.tibble(.)

### merge reads to clusters 
# merge reads to get reads regions (clusters)
reads_merged_gr <-
  reads_gr %>%
  unlist(.) %>%
  GenomicRanges::reduce(., ignore.strand = T)
reads_merged_gr$coordinates <- stringr::str_c(seqnames(reads_merged_gr), " ",
                                              start(reads_merged_gr), " ",
                                              end(reads_merged_gr))

# overlap reads with merged reads regions, get which read belongs to which cluster
reads_merged_df <- findOverlaps(reads_gr, reads_merged_gr, ignore.strand = TRUE) 
reads_merged_df <- tibble(read_id = reads_gr[queryHits(reads_merged_df)]$read_id, 
                          cluster_coord = reads_merged_gr[subjectHits(reads_merged_df)]$coordinates)

#### get genes info
genes_info <-
  gtf %>%
  as.data.frame(.) %>%
  as.tibble(.) %>%
  dplyr::distinct(., gene_id, .keep_all = T) %>%
  dplyr::select(gene_id, gene_name, gene_biotype)

#### overlap bam with exons, get width of overlap
reads_exons_overlaps <- findOverlaps(reads_gr, gtf, ignore.strand = TRUE)
reads_exons_filtered <- reads_gr[queryHits(reads_exons_overlaps)]
exons_filtered <- gtf[subjectHits(reads_exons_overlaps)]

# get width of overlap with exons, take maximum overlap
reads_exons_df <- 
  reads_exons_filtered %>% 
  as.data.frame(.) %>% 
  as.tibble(.) %>% 
  dplyr::mutate(exon_overlap = GenomicRanges::pintersect(reads_exons_filtered, exons_filtered, ignore.strand = TRUE) %>% width(.),
                gene_id = as.character(exons_filtered$gene_id)) %>% 
  dplyr::group_by(read_id) %>% 
  dplyr::filter(exon_overlap == max(exon_overlap)) %>%
  dplyr::ungroup(.) %>% 
  dplyr::distinct(read_id, .keep_all = T) %>% 
  dplyr::left_join(., genes_info, by = "gene_id") %>% 
  dplyr::select(read_id, exon_overlap, gene_id, gene_name, gene_biotype)
#### 

#### get width of total overlap between each read and repeatMasker
# find overlaps, filter ranges
reads_rmsk_gr_overlaps <- findOverlaps(reads_gr, rmsk_gr, ignore.strand = TRUE)
reads_rmsk_gr_filtered <- reads_gr[queryHits(reads_rmsk_gr_overlaps)]
rmsk_gr_filtered <- rmsk_gr[subjectHits(reads_rmsk_gr_overlaps)]

## get width of overlaps between reads and repeatMasker 
# (reduce because of overlaping features in repeatMasker)
rmsk_overlap_df <-
  GenomicRanges::pintersect(reads_rmsk_gr_filtered, rmsk_gr_filtered, ignore.strand = TRUE) %>% 
  GenomicRanges::split(., .$read_id) %>% 
  GenomicRanges::reduce(., ignore.strand = TRUE) %>% 
  GenomicRanges::width(.) %>% 
  sum(.) %>% 
  as.data.frame(.) %>% 
  tibble::rownames_to_column(., var = "read_id") %>% 
  as.tibble(.) %>% 
  data.table::setnames(., 2, "rmsk_overlap")
#### 

#### get repetitive elements which overlap start/end of read in sense direction
## start
# filter by overlap
reads_gr_start <- GenomicRanges::resize(reads_gr, width = 1, fix = "start", ignore.strand = T)
start_rmsk_overlaps <- findOverlaps(reads_gr_start, rmsk_gr, ignore.strand = TRUE)
reads_gr_start <- reads_gr_start[queryHits(start_rmsk_overlaps)]
rmsk_gr_start <- rmsk_gr[subjectHits(start_rmsk_overlaps)]

# get names/class/family of all repeats which overlap start
rmsk_start <- 
  tibble(read_id = reads_gr_start$read_id, 
         repName = rmsk_gr_start$repName, 
         repFamily = rmsk_gr_start$repFamily, 
         repClass = rmsk_gr_start$repClass, 
         rmsk_start_strand = as.character(strand(rmsk_gr_start))) %>% 
  dplyr::group_by(read_id) %>% 
  dplyr::summarise(repName = str_c(unique(repName), collapse = "|"), 
                   repFamily = str_c(unique(repFamily), collapse = "|"), 
                   repClass = str_c(unique(repClass), collapse = "|"), 
                   rmsk_start_strand = str_c(unique(rmsk_start_strand), collapse = "|")) %>% 
  dplyr::ungroup(.)

## end
# filter by overlap
reads_gr_end <- GenomicRanges::resize(reads_gr, width = 1, fix = "end", ignore.strand = T)
end_rmsk_overlaps <- findOverlaps(reads_gr_end, rmsk_gr, ignore.strand = TRUE)
reads_gr_end <- reads_gr_end[queryHits(end_rmsk_overlaps)]
rmsk_gr_end <- rmsk_gr[subjectHits(end_rmsk_overlaps)]

# get names/class/family of all repeats which overlap end
rmsk_end <- 
  tibble(read_id = reads_gr_end$read_id, 
         repName = rmsk_gr_end$repName, 
         repFamily = rmsk_gr_end$repFamily, 
         repClass = rmsk_gr_end$repClass, 
         rmsk_end_strand = as.character(strand(rmsk_gr_end))) %>% 
  dplyr::group_by(read_id) %>% 
  dplyr::summarise(repName = str_c(unique(repName), collapse = "|"), 
                   repFamily = str_c(unique(repFamily), collapse = "|"), 
                   repClass = str_c(unique(repClass), collapse = "|"), 
                   rmsk_end_strand = str_c(unique(rmsk_end_strand), collapse = "|")) %>% 
  dplyr::ungroup(.)

# find which reads are sense to repeatMasker elements 
# (repeats on + strand overlap with read start)
# (repeats on - strand overlap with read end) 
rmsk_sense_df <- 
  dplyr::full_join(rmsk_start %>% dplyr::select(read_id, rmsk_start_strand), 
                   rmsk_end %>% dplyr::select(read_id, rmsk_end_strand), 
                   by = "read_id") %>% 
  dplyr::mutate_at(.vars = c("rmsk_start_strand", "rmsk_end_strand"), .funs = funs(replace(., is.na(.), "no"))) %>% 
  dplyr::mutate(rmsk_sense = ifelse(((rmsk_start_strand == "+") & (rmsk_end_strand == "-")), "yes, ambiguous", "no")) %>% 
  dplyr::mutate(rmsk_sense = ifelse(((rmsk_start_strand == "+") & (rmsk_sense == "no")), "yes, plus strand", rmsk_sense)) %>% 
  dplyr::mutate(rmsk_sense = ifelse(((rmsk_end_strand == "-") & (rmsk_sense == "no")), "yes, minus strand", rmsk_sense)) %>% 
  dplyr::select(-c(rmsk_start_strand, rmsk_end_strand))

# filter start/end by sense 
rmsk_start %<>%
  dplyr::filter(read_id %in% (rmsk_sense_df %>% dplyr::filter((rmsk_sense == "yes, plus strand") | (rmsk_sense == "yes, ambiguous")) %$% read_id)) %>% 
  dplyr::select(-rmsk_start_strand)

rmsk_end %<>%
  dplyr::filter(read_id %in% (rmsk_sense_df %>% dplyr::filter((rmsk_sense == "yes, minus strand") | (rmsk_sense == "yes, ambiguous")) %$% read_id)) %>% 
  dplyr::select(-rmsk_end_strand)

# join together start/end features and sense data, join with total rmsk overlap
reads_rmsk_df <- 
  dplyr::full_join(rmsk_start, rmsk_end, by = "read_id") %>% 
  dplyr::mutate_at(.vars = vars(starts_with("rep")), .funs = funs(replace(., is.na(.), ""))) %>%
  dplyr::mutate(repName = str_c(repName.x, repName.y, sep = "||"), 
                repClass = str_c(repClass.x, repClass.y, sep = "||"), 
                repFamily = str_c(repFamily.x, repFamily.y, sep = "||")) %>% 
  dplyr::select(read_id, repName, repClass, repFamily) %>% 
  dplyr::full_join(., rmsk_sense_df, by = "read_id") %>% 
  dplyr::full_join(., rmsk_overlap_df, by = "read_id")


#### join all info about overlaping with exons and repeatMasker - clusters
class_clusters <-
  dplyr::full_join(reads_exons_df, reads_rmsk_df, by = "read_id") %>% 
  dplyr::full_join(., reads_df, by = "read_id") %>% 
  dplyr::left_join(., reads_merged_df, by = "read_id") %>%
  dplyr::mutate_at(.vars = c("exon_overlap", "rmsk_overlap"), .funs = funs(replace(., is.na(.), 0))) %>% 
  dplyr::mutate(exon_overlap = round((exon_overlap / width), digits = 3),
                rmsk_overlap = round((rmsk_overlap / width), digits = 3), 
                rmsk_sense = replace(rmsk_sense, is.na(rmsk_sense), "no")) %>%
  dplyr::mutate_at(.vars = vars(starts_with("rep")), .funs = funs(str_replace_all(., pattern = "\\|", ""))) %>% 
  dplyr::select(cluster_coord, exon_overlap_fraction = exon_overlap, 
                overlap_gene_name = gene_name, overlap_gene_biotype = gene_biotype,
                repeat_overlap_fraction = rmsk_overlap, sense_to_repeat = rmsk_sense, 
                sense_repName = repName, sense_repFamily = repFamily, sense_repClass = repClass) %>% 
  dplyr::arrange(cluster_coord) %T>% 
  readr::write_csv(x = ., path = file.path(outpath, "unspliced_reads_1C.clusters.csv"))

#### join all info about overlaping with exons and repeatMasker - individual reads
class_reads <-
  dplyr::full_join(reads_exons_df, reads_rmsk_df, by = "read_id") %>% 
  dplyr::full_join(., reads_df, by = "read_id") %>% 
  dplyr::mutate_at(.vars = c("exon_overlap", "rmsk_overlap"), .funs = funs(replace(., is.na(.), 0))) %>% 
  dplyr::mutate(exon_overlap = round((exon_overlap / width), digits = 3),
                rmsk_overlap = round((rmsk_overlap / width), digits = 3), 
                rmsk_sense = replace(rmsk_sense, is.na(rmsk_sense), "no")) %>%
  dplyr::mutate_at(.vars = vars(starts_with("rep")), .funs = funs(str_replace_all(., pattern = "\\|", ""))) %>% 
  dplyr::select(seqnames, start, end, read_length = width, 
                exon_overlap_fraction = exon_overlap, 
                overlap_gene_name = gene_name, overlap_gene_biotype = gene_biotype,
                repeat_overlap_fraction = rmsk_overlap, sense_to_repeat = rmsk_sense, 
                sense_repName = repName, sense_repFamily = repFamily, sense_repClass = repClass) %>% 
  dplyr::arrange(desc(read_length)) %T>% 
  readr::write_csv(x = ., path = file.path(outpath, "unspliced_reads_1C.width.csv"))
