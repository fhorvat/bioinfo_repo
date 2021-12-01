### INFO: counts small RNA reads mapped to different categories in repeatMasker and in ENSEMBL annotation (hierarchically)
### DATE: Fri May 25 08:37:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# mouse
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Analysis/class_reads/ES_DcrTrans_2012.SOLiD")
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
# source(file.path(lib_path, "wideScreen.R"))
# source(file.path(lib_path, "headt.R"))
# source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "bamToGRangesList.R"))
# wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# get gene info path
info_path <- list.files(genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.geneInfo.csv", full.names = T)

# get paths of reduced exons
exons_path <- list.files(genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.reducedExons.RDS", full.names = T)

# get repeatMasker path
rmsk_path <- list.files(genome_path, pattern = "rmsk.*[0-9]{6}.clean.fa.out.gz", full.names = T)

# get bam file path, name, experiment
bam_path <- "%BAM_PATH"
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Mapped/STAR_mm10.pCAG_EGFP_MosIR/s_RS7_NC_r1.SE.genome.Aligned.sortedByCoord.out.bam"
bam_name <- basename(bam_path) %>% str_remove_all(., ".SE.genome.Aligned.sortedByCoord.out.bam|.SE|.bam")
experiment_name <- "Eliska_mESC_MosIR"

######################################################## READ DATA
# read info about genes
genes_info <- readr::read_csv(info_path)

# gtf exons by genes
exons_gr <- readRDS(file = exons_path)

# read repeatMasker
rmsk_df <- readr::read_delim(rmsk_path, delim = "\t")

# read bam as ranges
bam_gr <- 
  bamToGRangesList(bam_path) %>% 
  unlist(.)

######################################################## MAIN CODE
### prepare ranges to count over - repeatMasker and exons of genes from ENSEMBL annotation
# get GRanges from repeatMasker
rmsk_gr <- 
  rmsk_df %>% 
  dplyr::mutate(gene_id = str_c(seqnames, ":", start, "-", end, "|", strand, "|", repName)) %>% 
  dplyr::select(seqnames:strand, gene_id, gene_biotype = repClass) %>% 
  dplyr::mutate(gene_biotype = replace(x = gene_biotype, 
                                       list = !(gene_biotype %in% c("rRNA", "LINE", "SINE", "LTR")), 
                                       values = "other_repeat")) %>% 
  GenomicRanges::GRanges(.)

# prepare exons
exons_gr_clean <-
  exons_gr %>% 
  as.data.frame(.) %>% 
  tibble::as.tibble(.) %>% 
  dplyr::left_join(., genes_info %>% dplyr::select(gene_id, gene_biotype), by = "gene_id") %>% 
  dplyr::select(seqnames:strand, gene_id, gene_biotype) %>% 
  GenomicRanges::GRanges(.)

# mos plasmid
mos_gr <- GenomicRanges::GRanges(seqnames = "pCAG-EGFP_MosIR", 
                                 ranges = IRanges(start = 1, end = 6660), 
                                 gene_id = "pCAG-EGFP_MosIR", 
                                 gene_biotype = "Mos_mRNA")

# join repeatMasker and ENSEMBL
count_ranges <- c(rmsk_gr, exons_gr_clean, mos_gr)


### class reads and save table
## hierarchy: Mos mRNA -> miRNA (sense only) -> mRNA (sense only) -> rRNA -> SINE -> LINE -> LTR -> other repeats -> annotated pseudogene -> other
# find overlaps between two GRanges - reads and repeatMasker/ENSEMBL/Mos plasmid
hits <- findOverlaps(bam_gr, count_ranges, ignore.strand = T)

# get hits in bam, convert to data.frame
read_hits <- 
  extractList(bam_gr, as(queryHits(hits), "List")) %>% 
  unlist(.)
read_hits$read_id <- names(read_hits)
names(read_hits) <- NULL
read_hits %<>% 
  as.data.frame(.) %>% 
  tibble::as.tibble(.) %>% 
  dplyr::mutate(read_coordinates = str_c(seqnames, start, end, sep = " ")) %>% 
  dplyr::select(read_id, read_coordinates, read_strand = strand)

# get hits in annotation, convert to data.frame
subject_hits <-
  extractList(count_ranges, as(subjectHits(hits), "List")) %>% 
  unlist(.) %>% 
  as.data.frame(.) %>% 
  as.tibble(.) %>% 
  dplyr::mutate(hit_coordinates = str_c(seqnames, start, end, sep = " ")) %>% 
  dplyr::select(gene_biotype, gene_id, hit_coordinates, hit_strand = strand)

# join tables
read_hits_all <- 
  dplyr::bind_cols(read_hits, subject_hits) %>% 
  dplyr::mutate(strand = ifelse(read_strand == hit_strand, "sense", "antisense"))

# get all biotypes to which one read maps
read_hits_sum <- 
  read_hits_all %>% 
  dplyr::mutate(gene_biotype = 
                  replace(gene_biotype, gene_biotype == "Mt_rRNA", "rRNA") %>% 
                  replace(., str_detect(., "pseudogene"), "annotated_pseudogene")) %>%
  dplyr::mutate(gene_biotype = ifelse(gene_biotype == "miRNA", str_c(gene_biotype, ".", strand), gene_biotype), 
                gene_biotype = ifelse(gene_biotype == "protein_coding", str_c(gene_biotype, ".", strand), gene_biotype)) %>% 
  dplyr::mutate(gene_biotype = replace(gene_biotype, 
                                       !(gene_biotype %in% c("Mos_mRNA", "miRNA.sense", "protein_coding.sense", "rRNA", "SINE", "LINE", "LTR", "other_repeat", "annotated_pseudogene")), 
                                       "other")) %>% 
  as.data.table(.)

# add read group
read_hits_sum <- read_hits_sum[ , read_group := ifelse("Mos_mRNA" %in% gene_biotype, "Mos_mRNA",
                                                       ifelse("miRNA.sense" %in% gene_biotype, "miRNA.sense",
                                                              ifelse("protein_coding.sense" %in% gene_biotype, "protein_coding.sense",
                                                                     ifelse("rRNA" %in% gene_biotype, "rRNA",
                                                                            ifelse("SINE" %in% gene_biotype, "SINE",
                                                                                   ifelse("LINE" %in% gene_biotype, "LINE",
                                                                                          ifelse("LTR" %in% gene_biotype, "LTR",
                                                                                                 ifelse("other_repeat" %in% gene_biotype, "other_repeat",
                                                                                                        ifelse("annotated_pseudogene" %in% gene_biotype, "annotated_pseudogene",
                                                                                                               ifelse("other" %in% gene_biotype, "other", "not_annotated")))))))))),
                                by = read_id]

# get and count reads not mapped to any category
reads_na <- 
  names(bam_gr) %>% 
  unique(.) %>% 
  .[!(. %in% read_hits_sum$read_id)] %>% 
  length(.)

# sum sense and antisense read groups
read_hits_final <- 
  read_hits_sum %>% 
  as.tibble(.) %>% 
  dplyr::select(read_id, read_group) %>% 
  unique(.) %>% 
  dplyr::group_by(read_group) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::bind_rows(., tibble(read_group = "not_annotated", count = reads_na)) %>% 
  dplyr::mutate(sample = bam_name, 
                experiment = experiment_name) %>% 
  dplyr::arrange(desc(count)) %T>% 
  readr::write_csv(x = ., path = file.path(outpath, str_c("smallRNA", experiment_name, bam_name, "csv", sep = ".")))
