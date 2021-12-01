### INFO: 
### DATE: 18. 07. 2017.  
### AUTHOR: Filip Horvat

# 1) fraction of the genome for each chromosome for mice
# 2) fraction of the mouse genome for 
# (i) protein-coding exons, 
# (ii) transcribed regions (exons+introns), 
# (iii) repetitive mobile elements total, 
# (iv) LINE1, 
# (v) LTR elements total, 
# (vi) ERVL elements

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/other")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(GenomicRanges)

######################################################## PATH VARIABLES
ensembl_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.gtf.gz"
# ensembl_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Ensembl_GRCm38.86.20161128.gtf.gz"

repeatmasker_path <- "/common/WORK/fhorvat/reference/mouse/mm10/UCSC/UCSC_repeatMaskerVIZ_OutBaseline_20160824.txt.gz"
# repeatmasker_path <- "/common/WORK/fhorvat/reference/mouse/mm10/UCSC/UCSC_repeatMasker_mm10_20170209_all_fields.txt.gz"

lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"

################################################################################### SOURCE FILES
source(file.path(lib_path, "GffToGRanges.R"))

######################################################## FUNCTIONS
getFraction <- function(genomic_ranges) {
  
  # takes GRanges, reduces and returns fraction of total genome size
  GenomicRanges::reduce(genomic_ranges, ignore.strand = T) %>%
    width(.) %>% 
    as.numeric(.) %>% 
    sum(.) %>% 
    magrittr::divide_by(genome_length) %>% 
    round(., digits = 5)
  
}

######################################################## READ DATA
### read in the ENSEMBL gene table
ensembl_gtf <- read_delim(file = ensembl_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c")) 

# read in repeatMasker table
rptmsk <- 
  read_delim(file = repeatmasker_path, delim = "\t") %>% 
  dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, repName, repClass, repFamily)

# read chromosome sizes
chr_size <-
  read_delim(file = "/common/WORK/fhorvat/reference/mouse/mm10/UCSC/mm10.chrom.sizes", delim = "\t", col_names = F) %>%
  magrittr::set_colnames(c("seqnames", "width")) %>% 
  dplyr::filter(!stringr::str_detect(string = seqnames, pattern = "Un|GL|JH"))

######################################################## MAIN CODE
### reformat gtf data.frame to GRanges
gtf_gr <- GffToGRanges(ensembl_gtf)

### get transcript with most exons for each gene
# get exons for each transcript
gtf_trans <- gtf_gr[gtf_gr$feature.type == "exon"]

# get unique gene_id/transcript_id combinations
gids <- unique(values(gtf_trans)[c('gene_id', 'transcript_id')])

# splits exon ranges on transcripts
gtf_trans <- split(gtf_trans, gtf_trans$transcript_id)

# orders transcripts based on number of exons in transcripts
gtf_trans <- gtf_trans[order(elementNROWS(gtf_trans), decreasing = T)]

# keeps only first transcript of each gene (the one with most exons)
gtf_trans <- gtf_trans[!duplicated(gids$gene_id[match(names(gtf_trans), gids$transcript_id)])]

### get total genome size
genome_length <- sum(as.numeric(chr_size$width))


##############
### 1) fraction of the genome for each chromosome for mice
chr_frac <- 
  chr_size %>% 
  dplyr::mutate(genome_fraction = round((width / genome_length), digits = 8)) %>% 
  dplyr::arrange(as.integer(stringr::str_replace(chr_frac$seqnames, "chr", ""))) %T>% 
  readr::write_csv(., path = "chromosome_frac.csv")

### 2) (i) protein-coding exons
protein_coding_exons <- 
  unlist(gtf_trans) %>%  
  .[(.$gene_biotype == "protein_coding") & !(duplicated(.$exon_id))] %>% 
  getFraction(.)

### 2) (ii) transcribed regions (exons + introns)
transcripts <- 
  gtf_gr[(gtf_gr$feature.type == "transcript") & (gtf_gr$transcript_id %in% names(gtf_trans))] %>% 
  getFraction(.)

### 2) (iii) repetitive mobile elements total
mobile_elements <- 
  rptmsk %>% 
  dplyr::filter(repClass %in% c("LINE", "SINE", "LTR")) %>% 
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T) %>% 
  getFraction(.)

### 2) (iv) LINE1
line1 <- 
  rptmsk %>% 
  dplyr::filter(repFamily == "L1") %>% 
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T) %>% 
  getFraction(.)

### 2)  (v) LTR elements total
ltr <- 
  rptmsk %>% 
  dplyr::filter(repClass == "LTR") %>% 
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T) %>% 
  getFraction(.)

### 2)  (vi) ERVL elements
ervl <- 
  rptmsk %>% 
  dplyr::filter(repFamily %in% c("ERVL", "ERVL-MaLR")) %>% 
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T) %>% 
  getFraction(.)

### put to data.frame, save
all_frac <- 
  rbind(protein_coding_exons, transcripts, mobile_elements, line1, ltr, ervl) %>% 
  as.data.frame(.) %>% 
  tibble::rownames_to_column("genomic_element") %>% 
  dplyr::rename(fraction = V1) %>% 
  readr::write_csv(., path = "genomic_elements_fraction_old.csv")