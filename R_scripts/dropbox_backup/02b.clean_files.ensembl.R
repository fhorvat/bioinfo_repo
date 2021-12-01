### INFO: creates table with relations between ensembl and USCS seqnames using assembly report from NCBI
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
# create vector of arguments from outside call
args <- commandArgs(TRUE)

# set ENSEMBL version
ensembl_release <- args[1]

# set ensembl species name
ensembl_name <- args[2]

# genome path
genome_path <- args[3]

# set working dir
setwd(genome_path)

### manual input
# ensembl_release <- 93
# ensembl_name <- "cchok1gshd"
# genome_path <- "/common/WORK/fhorvat/annotation/chinese_hamster/CHOK1GS_HDv1.GCA_900186095.1"
# setwd(genome_path)

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(biomaRt)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "gtfToGRanges.R"))

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# seqnames report path
seqnames_path <- list.files(path = inpath, pattern = ".*assembly_report.txt")

# ensembl .gtf path
ensembl_path <- list.files(path = inpath, pattern = str_c("ensembl.", ensembl_release, ".*[0-9]{6}.gtf.gz"))

# repeatMasker path
rmsk_path <- list.files(path = inpath, pattern = "rmsk.*raw.fa.out.gz")

######################################################## READ DATA
# read seqnames report
seqnames_table <- readr::read_delim(file = seqnames_path, delim = "\t", comment = "#", col_names = F)

# read gtf
ensembl_gtf <- read_delim(file = ensembl_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# read repeatMasker
rmsk_df <- readr::read_table2(file = rmsk_path, skip = 3, col_names = F)

######################################################## MAIN CODE
### clean and save repeatMasker
rmsk_df %>%
  dplyr::select(seqnames = X5, start = X6, end = X7, strand = X9, repName = X10, repClass_repFamily = X11) %>%
  tidyr::separate(col = repClass_repFamily, into = c("repClass", "repFamily"), sep = "/") %>%
  dplyr::mutate(strand = replace(strand, strand == "C", "-")) %T>%
  readr::write_delim(str_replace(rmsk_path, "raw", "clean"), delim = "\t")

# change gtf name
gtf_name <-
  basename(ensembl_path) %>% 
  stringr::str_remove(., ".gz")

### get additional info about genes from Biomart
## Ensembl versions
ensembl_url <-
  tibble(ens_version = c(93, 92, 91, 89, 86),
         date = c("Jul 2018", "Apr 2018", "Dec 2017", "May 2017", "Oct 2016"),
         URL_archive = c("http://jul2018.archive.ensembl.org", 
                         "http://apr2018.archive.ensembl.org",
                         "http://dec2017.archive.ensembl.org",
                         "http://may2017.archive.ensembl.org",
                         "http://oct2016.archive.ensembl.org")) %>%
  dplyr::filter(ens_version == ensembl_release) %$%
  URL_archive

# load Mart of mouse database from ensembl
# sometimes function useMart isn't able to connect to server and returns error, this chunck repeats useMart until there is no error
mart <- "error"
count <- 0
while(class(mart) == "character"){
  
  count <- count + 1
  print(str_c(ensembl_name, "_gene_ensembl", " ", count))
  
  # load ENSEMBL mart
  mart <- tryCatch(expr = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = stringr::str_c(ensembl_name, "_gene_ensembl"), host = ensembl_url),
                   error = function(x) return("error"))
  
  # if error try mouse strains mart
  if(class(mart) == "character"){
    mart <- tryCatch(expr = useMart(biomart = "ENSEMBL_MART_MOUSE", dataset = stringr::str_c(ensembl_name, "_gene_ensembl"), host = ensembl_url),
                     error = function(x) return("error"))
  }
  
  # stop if count get too big
  if(count > 2){
    stop("Something's not right")
  }
  
}

# get info about genes
ensembl_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "description"), mart = mart)

# merge with coordinates, save
gene_ID <-
  ensembl_gtf %>%
  gtfToGRanges(., filter = "gene") %>%
  as.data.frame(.) %>%
  as.tibble(.) %>%
  dplyr::select(gene_id, seqnames:end, strand) %>%
  dplyr::left_join(., ensembl_info, by = c("gene_id" = "ensembl_gene_id")) %>%
  dplyr::rename(gene_name = external_gene_name, gene_description = description) %T>%
  readr::write_csv(., path = file.path(outpath, gtf_name %>% stringr::str_replace(., ".gtf", ".geneInfo.csv")))

# get transcript - gene pairs, save
getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), mart = mart) %>%
  tibble::as.tibble(.) %>%
  dplyr::filter(ensembl_gene_id %in% gene_ID$gene_id) %>%
  dplyr::select(gene_id = ensembl_gene_id, transcript_id = ensembl_transcript_id) %>%
  dplyr::arrange(gene_id) %T>%
  readr::write_csv(., path = file.path(outpath, gtf_name %>% stringr::str_replace(., ".gtf", ".transcriptInfo.csv")))

### reduced exons
# convert GTF to GRanges, get only exons, reduce
exons_gr <-
  gtfToGRanges(ensembl_gtf, filter = "exon") %>%
  GenomicRanges::split(., .$gene_id) %>%
  GenomicRanges::reduce(., ignore.strand = T) %>%
  unlist(.)

# add strand info, split again
exons_gr$gene_id <- names(exons_gr)
names(exons_gr) <- NULL
exons_gr <-
  exons_gr %>%
  as.data.frame(.) %>%
  as.tibble(.) %>%
  dplyr::select(-strand) %>%
  dplyr::left_join(., gene_ID %>% dplyr::select(gene_id, strand), by = "gene_id") %>%
  GenomicRanges::GRanges(.) %>%
  GenomicRanges::split(., .$gene_id)

# save RDS
saveRDS(object = exons_gr, file = file.path(outpath, gtf_name %>% stringr::str_replace(., ".gtf", ".reducedExons.RDS")))
