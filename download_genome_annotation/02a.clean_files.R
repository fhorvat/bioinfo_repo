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
#ensembl_release <- 96
#ensembl_name <- "btaurus"
#genome_path <- "/common/DB/genome_reference/cow/bosTau9.ARS-UCD1.2.GCA_002263795.2"
#setwd(genome_path)

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

# UCSC fasta index path
UCSC_faidx_path <- list.files(path = inpath, pattern = "\\.fai")

######################################################## READ DATA
# read seqnames report
seqnames_table <- readr::read_delim(file = seqnames_path, delim = "\t", comment = "#", col_names = F)

# read gtf
ensembl_gtf <- read_delim(file = ensembl_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# read repeatMasker
rmsk_df <- readr::read_table2(file = rmsk_path, skip = 3, col_names = F)

# read UCSC fasta index
UCSC_faidx <-
  readr::read_delim(UCSC_faidx_path, delim = "\t", col_names = F) %>%
  dplyr::select(UCSC_name = X1)

######################################################## MAIN CODE
### clean and save repeatMasker
rmsk_df %>%
  dplyr::select(seqnames = X5, start = X6, end = X7, strand = X9, repName = X10, repClass_repFamily = X11, rmsk_id = X15) %>%
  tidyr::separate(col = repClass_repFamily, into = c("repClass", "repFamily"), sep = "/") %>%
  dplyr::mutate(strand = replace(strand, strand == "C", "-")) %T>%
  readr::write_delim(str_replace(rmsk_path, "raw", "clean"), delim = "\t")

### UCSC seqnames
# For unplaced and unlocalized scaffolds  UCSC almost always uses genBank accession = genBank_accn (for example mouse),
# but sometimes it uses refSeq accession = refSeq_accn (for example pig, cow).
# Sometimes table already has UCSC-style names.
# ENSEMBL uses GenBank accession for unplaced and unlocalized scaffolds.
UCSC_seqnames <-
  seqnames_table %>%
  magrittr::set_colnames(., c("sequence_name", "sequence_role", "assigned_molecule",
                              "assigned_molecule_location_type", "genBank_accn",
                              "relationship", "refSeq_accn", "assembly_unit",
                              "sequence_length", "UCSC_style_name")) %>%
  dplyr::select(sequence_role, assigned_molecule, genBank_accn, sequence_name, assigned_molecule_location_type, refSeq_accn, UCSC_style_name) %>%
  dplyr::mutate(ensembl_name = ifelse(sequence_role == "assembled-molecule",
                                      assigned_molecule,
                                      genBank_accn)) %>%
  dplyr::mutate(UCSC_name = ifelse(sequence_role == "assembled-molecule",
                                   str_c("chr", assigned_molecule),
                                   str_c("chr", assigned_molecule, "_", genBank_accn %>% str_replace(., "\\.", "v"))),
                UCSC_name = ifelse(sequence_role == "unlocalized-scaffold", str_c(UCSC_name, "_random"), UCSC_name),
                UCSC_name = replace(UCSC_name, assigned_molecule_location_type == "Mitochondrion", "chrM"),
                UCSC_name = str_replace(UCSC_name, "chrna", "chrUn")) %>%
  dplyr::mutate(UCSC_name = ifelse(UCSC_style_name == "na", UCSC_name, UCSC_style_name)) %>%
  dplyr::select(ensembl_name, UCSC_name, sequence_role) %T>%
  readr::write_delim(x = .,
                     path = file.path(outpath, seqnames_path %>% stringr::str_replace(., "_assembly_report.txt", ".ensembl2UCSC.txt")),
                     delim = "\t")

# sanity check
ensembl_names <- all(unique(ensembl_gtf$X1) %in% UCSC_seqnames$ensembl_name)
UCSC_names <- all(UCSC_faidx$UCSC_name %in% UCSC_seqnames$UCSC_name)
if(!(ensembl_names & UCSC_names)){
  stop("Ensembl: ", ensembl_names, "; UCSC: ", UCSC_names)
}

### ensembl to UCSC .gtf
# decide whether to remove or keep scaffolds from assembly
primary <- F

# change gtf name
gtf_name <-
  basename(ensembl_path) %>%
  stringr::str_replace(., ".gtf.gz", ".UCSCseqnames.gtf") %>%
  {if(primary) stringr::str_replace(., ".gtf", ".primary.gtf") else .}

# join gtf with UCSC seqnames, optional keep only primary assembly, write as .gtf
ensembl_gtf_UCSC <-
  ensembl_gtf %>%
  dplyr::left_join(., UCSC_seqnames, by = c("X1" = "ensembl_name")) %>%
  {if(primary) dplyr::filter(., sequence_role == "assembled-molecule") else .} %>%
  dplyr::select(X1 = UCSC_name, X2:X9) %T>%
  write.table(x = ., file = file.path(outpath, gtf_name), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

# gzip gtf
system(stringr::str_c("gzip ", file.path(outpath, gtf_name)))


### get additional info about genes from Biomart
## Ensembl versions
ensembl_url <-
  tibble(ens_version = c(99, 98, 96, 95, 94, 93, 92, 91, 89, 86),
         date = c("Jan 2020", "Sep 2019", "Apr 2019", "Jan 2019", "Oct 2018", "Jul 2018", "Apr 2018", "Dec 2017", "May 2017", "Oct 2016"),
         URL_archive = c("www.ensembl.org",
                         "http://sep2019.archive.ensembl.org",
                         "http://apr2019.archive.ensembl.org",
                         "http://jan2019.archive.ensembl.org",
                         "http://oct2018.archive.ensembl.org",
                         "http://jul2018.archive.ensembl.org",
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
  ensembl_gtf_UCSC %>%
  gtfToGRanges(., filter = "gene") %>%
  as.data.frame(.) %>%
  as_tibble(.) %>%
  dplyr::select(gene_id, seqnames:end, strand) %>%
  dplyr::left_join(., ensembl_info, by = c("gene_id" = "ensembl_gene_id")) %>%
  dplyr::rename(gene_name = external_gene_name, gene_description = description) %T>%
  readr::write_csv(., path = file.path(outpath, gtf_name %>% stringr::str_replace(., ".gtf", ".geneInfo.csv")))

# get transcript - gene pairs, save
getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), mart = mart) %>%
  tibble::as_tibble(.) %>%
  dplyr::filter(ensembl_gene_id %in% gene_ID$gene_id) %>%
  dplyr::select(gene_id = ensembl_gene_id, transcript_id = ensembl_transcript_id) %>%
  dplyr::arrange(gene_id) %T>%
  readr::write_csv(., path = file.path(outpath, gtf_name %>% stringr::str_replace(., ".gtf", ".transcriptInfo.csv")))

### reduced exons
# convert GTF to GRanges, get only exons, reduce
exons_gr <-
  gtfToGRanges(ensembl_gtf_UCSC, filter = "exon") %>%
  GenomicRanges::split(., .$gene_id) %>%
  GenomicRanges::reduce(., ignore.strand = T) %>%
  unlist(.)

# add strand info, split again
exons_gr$gene_id <- names(exons_gr)
names(exons_gr) <- NULL
exons_gr <-
  exons_gr %>%
  as.data.frame(.) %>%
  as_tibble(.) %>%
  dplyr::select(-strand) %>%
  dplyr::left_join(., gene_ID %>% dplyr::select(gene_id, strand), by = "gene_id") %>%
  GenomicRanges::GRanges(.) %>%
  GenomicRanges::split(., .$gene_id)

# save RDS
saveRDS(object = exons_gr, file = file.path(outpath, gtf_name %>% stringr::str_replace(., ".gtf", ".reducedExons.RDS")))
