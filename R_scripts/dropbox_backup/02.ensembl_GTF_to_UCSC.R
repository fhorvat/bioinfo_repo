### INFO: creates table with relations between ensembl and USCS seqnames using assembly report from NCBI
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
# set ENSEMBL version
ens_ver <- 91

# set ensembl species name
ensembl_species <- "mmusculus"

# genome path
genome_path <- "/common/DB/genome_reference/cow/bosTau8.UMD3.1.GCA_000003055.4"

# set working dir
setwd(genome_path)

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
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "gtfToGRanges.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

#seqnames report path
seqnames_path <- list.files(path = inpath, pattern = ".*assembly_report.txt")

# ensembl .gtf path
ensembl_path <- list.files(path = inpath, pattern = str_c("ensembl.", ens_ver, ".*[0-9]{6}.gtf.gz"))

######################################################## READ DATA
# read seqnames report
seqnames_table <- readr::read_delim(file = seqnames_path, delim = "\t", comment = "#", col_names = F)

# read gtf
ensembl_gtf <- read_delim(file = ensembl_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

######################################################## MAIN CODE
### UCSC seqnames 
# For unplaced and unlocalized scaffolds  UCSC almost always uses genBank accession = genBank_accn (for example cow, mouse), 
# but sometimes it uses refSeq accession = refSeq_accn (for example pig). 
# Sometimes table already has UCSC-style names. 
# ENSEMBL uses GenBank accession for unplaced and unlocalized scaffolds. 
USCS_seqnames <-
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
  tibble(ens_version = c(92, 91, 86), 
         date = c("Apr 2018", "91 Dec 2017", "Oct 2016"), 
         URL_archive = c("http://apr2018.archive.ensembl.org", 
                         "http://dec2017.archive.ensembl.org", 
                         "http://oct2016.archive.ensembl.org")) %>% 
  dplyr::filter(ens_version == ens_ver) %$% 
  URL_archive

# load Mart of mouse database from ensembl
# sometimes function useMart isn't able to connect to server and returns error, this chunck repeats useMart until there is no error
mart <- "error"
count <- 0
while(class(mart) == "character"){
  count <- count + 1
  print(str_c(ensembl_species, "_gene_ensembl", " ", count))
  mart <- tryCatch(expr = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = str_c(ensembl_species, "_gene_ensembl"), 
                                  host = "http://dec2017.archive.ensembl.org"), 
                   error = function(x) return("error"))
}

# get info about genes, merge
ensembl_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "description"), mart = mart)

# get gene ID's and other info from biomart, save
gene_ID <- 
  ensembl_gtf_UCSC %>% 
  gtfToGRanges(., filter = "gene") %>% 
  as.data.frame(.) %>% 
  as.tibble(.) %>% 
  dplyr::select(gene_id, seqnames:end, strand) %>% 
  dplyr::left_join(., ensembl_info, by = c("gene_id" = "ensembl_gene_id")) %>% 
  dplyr::rename(gene_name = external_gene_name, gene_description = description) %T>%
  readr::write_csv(., path = file.path(outpath, gtf_name %>% stringr::str_replace(., ".gtf", ".geneInfo.csv")))


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
  as.tibble(.) %>% 
  dplyr::select(-strand) %>% 
  dplyr::left_join(., gene_ID %>% dplyr::select(gene_id, strand), by = "gene_id") %>% 
  GenomicRanges::GRanges(.) %>% 
  GenomicRanges::split(., .$gene_id)

# save RDS
saveRDS(object = exons_gr, file = file.path(outpath, gtf_name %>% stringr::str_replace(., ".gtf", ".reducedExons.RDS")))
