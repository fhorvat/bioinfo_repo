### INFO: counts small RNA reads mapped to different categories in repeatMasker and in ENSEMBL annotation (hierarchically)
### DATE: Fri May 25 08:37:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# mouse
# setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/mouse")
# genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# cow
# setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/cow")
# genome_path <- "/common/DB/genome_reference/cow/bosTau8.UMD3.1.GCA_000003055.4"

# pig
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/pig")
genome_path <- "/common/DB/genome_reference/pig/susScr11.Sscrofa11.1.GCA_000003025.6"

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "bamToGRangesList.R"))

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
bam_name <- basename(bam_path) %>% str_remove(., ".SE.perfect.bam")
experiment_name <- bam_path %>% str_remove(., "/Filtered/.*") %>% basename(.) %>% str_remove(., "(?<=[0-9])_.*")

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
# get reads with width between 21 and 23 nt
bam_gr_filt <- bam_gr[(width(bam_gr) >= 21) & (width(bam_gr) <= 23)]


## class reads - rRNA -> other -> LINE -> SINE -> LTR -> miRNA -> lincRNA ->  protein_coding
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
  dplyr::mutate(gene_biotype = replace(x = gene_biotype,
                                       list = gene_biotype == "Mt_rRNA" ,
                                       values = "rRNA"),
                gene_biotype = replace(x = gene_biotype,
                                       list = str_detect(gene_biotype, "pseudogene"),
                                       values = "annotated_pseudogene"),
                gene_biotype = replace(x = gene_biotype,
                                       list = !(gene_biotype %in% c("rRNA", "SINE", "LINE", "LTR", "other_repeat", "miRNA", "lincRNA", "protein_coding", "annotated_pseudogene")),
                                       values = "other")) %>%
  dplyr::group_by(read_id, strand) %>%
  dplyr::summarise(gene_biotype = str_c(unique(gene_biotype), collapse = ", "))

# divide to sense and antisense
read_hits_sense <-
  read_hits_sum %>%
  dplyr::filter(strand == "sense")

read_hits_antisense <-
  read_hits_sum %>%
  dplyr::filter(strand == "antisense")

# hierarchically categorize reads in loop
reads_sense <- tibble()
reads_antisense <- tibble()

for(biotype_class in c("rRNA", "LINE", "SINE", "LTR", "other_repeat", "miRNA", "lincRNA", "protein_coding", "annotated_pseudogene", "other")){
  
  ### sense
  # filter data.frame
  class_reads <-
    read_hits_sense %>%
    dplyr::filter(str_detect(string = gene_biotype, pattern = biotype_class)) %>%
    dplyr::mutate(gene_biotype = biotype_class)
  
  # add column to sum table
  reads_sense %<>%
    dplyr::bind_rows(., class_reads)
  
  # filter original table
  read_hits_sense %<>%
    dplyr::filter(!(read_id %in% class_reads$read_id))
  
  
  ### antisense
  # filter data.frame
  class_reads <-
    read_hits_antisense %>%
    dplyr::filter(str_detect(string = gene_biotype, pattern = biotype_class)) %>%
    dplyr::mutate(gene_biotype = biotype_class)
  
  # add column to sum table
  reads_antisense %<>%
    dplyr::bind_rows(., class_reads)
  
  # filter original table
  read_hits_antisense %<>%
    dplyr::filter(!(read_id %in% class_reads$read_id))
  
}

# get and count reads not mapped to any category
reads_na <-
  names(bam_gr_filt) %>%
  unique(.) %>%
  .[!(. %in% c(reads_sense$read_id, reads_antisense$read_id))] %>%
  unique(.) %>%
  length(.)

# join reads category with strand data, sum sense and antisense
read_hits_final <-
  rbind(reads_sense, reads_antisense) %>%
  dplyr::group_by(gene_biotype, strand) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::bind_rows(., tibble(gene_biotype = "not_annotated", strand = NA, count = reads_na)) %>%
  dplyr::mutate(sample = bam_name,
                experiment = experiment_name) %T>%
  readr::write_csv(x = ., path = file.path(outpath, str_c("smallRNA.perfect.21to23.summary.counts", experiment_name, bam_name, "csv", sep = ".")))

