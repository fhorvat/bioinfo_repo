### INFO: creates table with relations between ensembl and USCS seqnames using assembly report from NCBI
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
## manual input
genome_path <- "/common/DB/genome_reference/rat/rn7.mRatBN7.2.GCA_015227675.2"
setwd(genome_path)

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# seqnames report path
seqnames_path <- list.files(path = inpath, pattern = ".*assembly_report.txt")

# repeatMasker path
rmsk_path <- list.files(path = inpath, pattern = "rmsk.*raw.fa.out.gz")

# refSeq .gtf path
gtf_path <- list.files(path = inpath, pattern = ".*\\.gff|.*\\.gff.gz", full.names = T)
gtf_path <- gtf_path[!str_detect(gtf_path, "UCSC")]

# UCSC fasta index path
UCSC_faidx_path <- list.files(path = inpath, pattern = "\\.fai")

######################################################## READ DATA
# read seqnames report
seqnames_table <- readr::read_delim(file = seqnames_path, delim = "\t", comment = "#", col_names = F)

# read repeatMasker
rmsk_df <- readr::read_table2(file = rmsk_path, skip = 3, col_names = F)

# read refSeq .gtf as tibble
gtf_tb <- read_delim(file = gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# read UCSC fasta index
UCSC_faidx <-
  readr::read_delim(UCSC_faidx_path, delim = "\t", col_names = F) %>%
  dplyr::select(UCSC_name = X1)

######################################################## MAIN CODE
### clean and save repeatMasker
rmsk_df_tidy <- 
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
  dplyr::select(sequence_role, assigned_molecule, genBank_accn, sequence_name, 
                assigned_molecule_location_type, refSeq_accn, UCSC_style_name) %>%
  dplyr::mutate(UCSC_name = ifelse(sequence_role == "assembled-molecule",
                                   str_c("chr", assigned_molecule),
                                   str_c("chr", assigned_molecule, "_", genBank_accn %>% str_replace(., "\\.", "v"))),
                UCSC_name = ifelse(sequence_role == "unlocalized-scaffold", str_c(UCSC_name, "_random"), UCSC_name),
                UCSC_name = replace(UCSC_name, assigned_molecule_location_type == "Mitochondrion", "chrM"),
                UCSC_name = str_replace(UCSC_name, "chrna", "chrUn")) %>%
  dplyr::mutate(UCSC_name = ifelse(UCSC_style_name == "na", UCSC_name, UCSC_style_name)) %>%
  dplyr::select(refSeq_accn, UCSC_name, sequence_role) %T>%
  readr::write_delim(x = .,
                     file = file.path(outpath, seqnames_path %>% stringr::str_replace(., "_assembly_report.txt", ".refSeq2UCSC.txt")),
                     delim = "\t")

# change gtf name
gtf_name <-
  basename(gtf_path) %>%
  stringr::str_replace(., ".gff.gz", ".UCSCseqnames.gff")

# join gtf with UCSC seqnames
gtf_tb_UCSC <-
  gtf_tb %>%
  dplyr::left_join(., UCSC_seqnames %>% dplyr::select(-sequence_role), by = c("X1" = "refSeq_accn")) %>%
  dplyr::select(X1 = UCSC_name, X2:X9)

# check
if(gtf_name == basename(gtf_path)){
  
  # stop
  stop("You're gonna overwrite your original .gff!")
  
}else{
  
  # write 
  write.table(x = gtf_tb_UCSC, file = file.path(outpath, gtf_name), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  # gzip gtf
  system(stringr::str_c("gzip ", file.path(outpath, gtf_name)))
  
}




### get gene info
# refSeq .gtf path
gtf_path <- list.files(path = inpath, pattern = ".*\\.UCSCseqnames.gff|.*\\.UCSCseqnames.gff.gz", full.names = T)

# read refSeq .gtf as GRanges
gtf_gr <- rtracklayer::import(gtf_path)

# get genes
gtf_gr_gene <- gtf_gr[mcols(gtf_gr)$type %in% c("gene", "pseudogene")]

# create tibble
refseq_geneInfo <- tibble(gene_id = mcols(gtf_gr_gene)$gene, 
                          seqnames = as.character(seqnames(gtf_gr_gene)),
                          start = as.numeric(start(gtf_gr_gene)), 
                          end = as.numeric(end(gtf_gr_gene)), 
                          strand = as.character(strand(gtf_gr_gene)), 
                          gene_name =  as.character(mcols(gtf_gr_gene)$gene), 
                          gene_biotype = as.character(mcols(gtf_gr_gene)$gene_biotype), 
                          gene_description = mcols(gtf_gr_gene)$description)


### reduced exons
# get clean table
gtf_tb_tidy <-
  gtf_gr %>% 
  as_tibble(.) %>% 
  dplyr::select(seqnames:strand, type, ID, Name, gene, gene_biotype, Parent)

# get parents
parents <- 
  gtf_tb_tidy$Parent %>% 
  purrr::map(., function(x) replace(x, length(x) == 0, NA)) %>% 
  unlist(.)

# add to table
exons_tb <- 
  gtf_tb_tidy %>% 
  dplyr::mutate(Parent = parents) %>% 
  dplyr::filter(type == "exon")

# add mitochondrion genes to table
exons_mt <- 
  gtf_tb_tidy %>% 
  dplyr::filter(seqnames == "chrM", type == "gene")

# add some pseudogenes which don't have annotated exons to table
pseudogenes_list <- 
  refseq_geneInfo %>% 
  dplyr::filter(!(gene_id %in% exons_tb$gene), 
                gene_biotype == "pseudogene") %$%
  gene_name %>% 
  unique(.)

# filter
exons_ps <-
  gtf_tb_tidy %>% 
  dplyr::filter(gene %in% pseudogenes_list)

# convert GTF to GRanges, get only exons, reduce
exons_gr <-
  rbind(exons_tb, exons_mt, exons_ps) %>% 
  GRanges(.) %>% 
  GenomicRanges::split(., .$gene) %>%
  GenomicRanges::reduce(., ignore.strand = T) %>%
  unlist(.)

# add names as gene ID
exons_gr$gene_id <- names(exons_gr)
names(exons_gr) <- NULL

# add strand info, split again
exons_gr_strand <-
  exons_gr %>%
  as.data.frame(.) %>%
  as_tibble(.) %>%
  dplyr::select(-strand) %>%
  dplyr::left_join(., refseq_geneInfo %>% dplyr::select(gene_id, strand), by = "gene_id") %>%
  GenomicRanges::GRanges(.) %>%
  GenomicRanges::split(., .$gene_id)

# get name
out_name <- 
  basename(gtf_path) %>% 
  str_replace(., "\\.gff$|\\.gff.gz$", ".reducedExons.RDS")

# check
if(out_name == basename(gtf_path)){
  
  # stop
  stop("Don't overwrite your .gff file!")
  
}else{
  
  # save RDS
  saveRDS(object = exons_gr_strand, file = file.path(outpath, out_name))
  
}
