### INFO: counts small RNA reads mapped to different categories in repeatMasker and in ENSEMBL annotation (hierarchically)
### DATE: Fri Jul 20 23:46:48 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# mouse
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Analysis/class_reads/detailed_other")
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

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "bamToGRangesList.R"))
wideScreen()

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
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Mapped/STAR_mm10_new/s_RSP_Mos_r3.SE.genome.Aligned.sortedByCoord.out.bam"
bam_name <- basename(bam_path) %>% str_remove(., ".SE.genome.Aligned.sortedByCoord.out.bam")
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

# join repeatMasker and ENSEMBL
count_ranges <- c(rmsk_gr, exons_gr_clean)


### prepare bam ranges to count
# get all reads
bam_gr_filt <- bam_gr[(width(bam_gr) >= 21) & (width(bam_gr) <= 23)]


### class reads
## hierarchy: Mos mRNA -> miRNA (sense only) -> mRNA (sense only) -> rRNA -> SINE -> LINE -> LTR -> other repeats -> annotated pseudogene -> other
# find overlaps between two GRanges - reads and repeatMasker/ENSEMBL
hits <- findOverlaps(bam_gr_filt, count_ranges, ignore.strand = T)

# get hits in bam, convert to data.frame
read_hits <-
  extractList(bam_gr_filt, as(queryHits(hits), "List")) %>%
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
                  replace(., gene_id == "ENSMUSG00000078365", "Mos_mRNA") %>%
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


##### test other categories
# get read IDs which are in "other_repeat" and "other" category
reads_filtered <- 
  read_hits_sum %>%
  as.tibble(.) %>%
  dplyr::select(read_id, read_group) %>%
  unique(.) %>% 
  dplyr::filter(read_group == "other" | read_group == "other_repeat") %$% 
  read_id

# get sequences of reads in those categories
bam_other_filt <- 
  GenomicAlignments::readGAlignmentsList(bam_path, use.names = TRUE, param = ScanBamParam(what = c("qname", "seq"))) %>% 
  unlist(.) %>% 
  .[names(.) %in% reads_filtered]

# get sequences
bam_seq <- 
  mcols(bam_other_filt)$seq %>% 
  as.character(.) %>% 
  tibble(seq = .) %>% 
  dplyr::group_by(seq) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::top_n(n = 30, wt = count) %>% 
  dplyr::arrange(desc(count))

### overlap with full repeatMasker
# get GRanges from repeatMasker
rmsk_gr <-
  rmsk_df %>%
  dplyr::mutate(gene_id = str_c(seqnames, ":", start, "-", end, "|", strand, "|", repName)) %>%
  dplyr::select(seqnames:strand, gene_id, gene_biotype = repClass) %>% 
  GenomicRanges::GRanges(.)

# join repeatMasker and ENSEMBL
count_ranges <- c(rmsk_gr, exons_gr_clean)

# find overlaps between two GRanges - reads and repeatMasker/ENSEMBL
hits <- findOverlaps(bam_other, count_ranges, ignore.strand = T)

# get hits in bam, convert to data.frame
read_hits <-
  extractList(bam_other, as(queryHits(hits), "List")) %>%
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

# join tables, sum gene biotypes
read_hits_all <-
  dplyr::bind_cols(read_hits, subject_hits) %>%
  dplyr::mutate(strand = ifelse(read_strand == hit_strand, "sense", "antisense")) %>% 
  dplyr::group_by(gene_biotype) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::arrange(desc(count)) %>% 
  dplyr::mutate(sample = bam_name, 
                experiment = experiment_name) %T>% 
  readr::write_csv(x = ., path = file.path(outpath, str_c("smallRNA.all.21to23", experiment_name, bam_name, "detailed_other.csv", sep = ".")))

# # # join tables, sum gene ID's
# read_hits_all <-
#   dplyr::bind_cols(read_hits, subject_hits) %>%
#   dplyr::mutate(strand = ifelse(read_strand == hit_strand, "sense", "antisense")) %>% 
#   dplyr::group_by(gene_id) %>% 
#   dplyr::summarise(count = n()) %>% 
#   dplyr::arrange(desc(count)) %>% 
#   dplyr::left_join(genes_info, by = "gene_id") %>% 
#   tidyr::unite(col = "coordinates", seqnames, start, end, strand, sep = " ") %T>% 
#   readr::write_csv(x = ., path = file.path(outpath, str_c("smallRNA.all_widths", experiment_name, bam_name, "detailed_other.genes.csv", sep = ".")))
