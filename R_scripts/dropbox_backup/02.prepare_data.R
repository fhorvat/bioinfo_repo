### INFO: 
### DATE: Sun Nov 04 19:08:17 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/accessory_data")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)
library(biomaRt)
library(targetscan.Mm.eg.db)
library(seqinr)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "gtfToGRanges.R"))

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# # genome path
# genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# # ensembl .gtf path
# ensembl_path <- list.files(path = genome_path, pattern = "ensembl.93.*[0-9]{6}.UCSC.*gtf.gz", full.names = T)

# # multiple alignment of 3' UTRs between mouse and rat/dog/human
# aligned_UTRs_path <- list.files(file.path(inpath, "multiz60way.3primeUTRs"), ".*3primeUTRs.txt", full.names = T)

# all 3'UTR sequences 
all_3UTRs_path <- file.path(inpath, "targetScan", "targetScan.mouse.7.1.20181112.UTR_sequences.txt.gz")

# taxonomy ID path
tax_id_path <- file.path(inpath, "ncbi_taxonomy_names.clean.txt")

######################################################## READ DATA
# # read ensembl .gtf
# ensembl_gtf <- read_delim(file = ensembl_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# # read aligned UTRs
# aligned_UTRs <- purrr::map(aligned_UTRs_path, function(path){
#   
#   # get target name
#   target_animal <- str_extract(path, "rn5|canFam3|hg19")
#   
#   # fread table
#   aligned_dt <- fread(path, sep = " ", col.names = c("mouse_coords_seq", str_c(target_animal, "_coords"), str_c(target_animal, "_seq")))
#   setkey(aligned_dt, "mouse_coords_seq")
#   
# }) %>% 
#   purrr::reduce(., merge, all = T) %>% 
#   .[complete.cases(.)]

# read all aligned 3'UTRs from targetScan
all_3UTRs <-  fread(all_3UTRs_path, col.names = c("refseq_id", "gene_id", "gene_symbol", "tax_id", "UTR_seq"), colClasses = "character")
all_3UTRs[tax_id == "9478", tax_id := "1868482"]
setkey(all_3UTRs, tax_id)

# read tax ID
tax_id <- fread(tax_id_path, col.names = c("tax_id", "tax_name", "unique_name", "name_class"), quote = "", colClasses = "character")
setkey(tax_id, tax_id)

######################################################## MAIN CODE
#### DISTANCE GRID ####
# prepare full grid, save
purrr::map(c(61, 21, 50), function(dimension){
  
  as.tibble(expand.grid(1:dimension, 1:dimension)) %>% 
    set_colnames(c("bin_absolute", "bin_relative")) %>% 
    dplyr::mutate(bin_id = str_c(bin_absolute, bin_relative, sep = ".")) %>% 
    column_to_rownames("bin_id") %>%
    dist(., method = "euclidean") %>% 
    as.matrix(.) %>%
    as.data.table(., keep.rownames = "bin_id") %>% 
    melt(., id.vars = "bin_id", variable.name = "bin_id_dist", value.name = "dist") %T>%
    saveRDS(., file = file.path(outpath, str_c("distance_grid.Euclidean.", dimension, "x", dimension, ".dt.RDS")))
  
  # return
  return(dimension)
  
})


#### SEPTAMERS ####
## expand 7mers
kmer_7 <- 
  expand.grid(rep(list(c("A", "G", "T", "C")), 7)) %>% 
  do.call(str_c, .)

# create table, save
kmer_7_tb <- 
  tibble(kmer = kmer_7) %>% 
  mutate(seed_7mer_m8 = DNAStringSet(kmer_7) %>% reverseComplement(.) %>% as.character(.), 
         seed_7mer_1a = str_c(str_sub(seed_7mer_m8, 2, 7), "A")) %>% 
  dplyr::rowwise(.) %>% 
  mutate_at(vars(contains("seed")), funs(pattern = str_c(str_split(., "", simplify = T), collapse = "_*"))) %T>%
  saveRDS(., file.path(outpath, "seeds.7mers_patterns.RDS"))


#### ALIGNED 3' UTRs ####
## filter aligned 3'UTRs from targetScan (get only mouse, rat, cow and human) and save as fasta and in table
# clean tax_id table
tax_id %<>%
  .[tax_id %in% unique(all_3UTRs$tax_id) & name_class == "scientific name", .(tax_id, tax_name)] %>% 
  .[, tax_name := str_replace_all(tax_name, " ", "_")] %>% 
  .[]

# add taxon name to sequence table
all_3UTRs[tax_id, tax_name := i.tax_name, on = "tax_id"]

# get only human, rat, mouse and cow UTRs
filtered_3UTRs <- all_3UTRs[tax_name %in% c("Bos_taurus", "Homo_sapiens", "Mus_musculus", "Rattus_norvegicus")] 

# add tax name to gene_id, replace U with T
filtered_3UTRs[, `:=`(gene_id = str_c(gene_id, "|", tax_name), 
                      UTR_seq = str_replace_all(UTR_seq, "U|u", "T"))]

# convert to DNAStringSet, save
filtered_3UTRs <-
  DNAStringSet(filtered_3UTRs$UTR_seq) %>% 
  set_names(., filtered_3UTRs$gene_id) %T>% 
  writeXStringSet(., filepath = file.path(outpath, "targetScan", "targetScan.mouse.7.1.20181112.UTR_sequences.filtered.fasta"))

# also save as list of DNAStringSets
filtered_3UTRs_list <- 
  list("Mus_musculus" = filtered_3UTRs[str_detect(names(filtered_3UTRs), "Mus_musculus")], 
       "Rattus_norvegicus" = filtered_3UTRs[str_detect(names(filtered_3UTRs), "Rattus_norvegicus")], 
       "Bos_taurus" = filtered_3UTRs[str_detect(names(filtered_3UTRs), "Bos_taurus")], 
       "Homo_sapiens" = filtered_3UTRs[str_detect(names(filtered_3UTRs), "Homo_sapiens")]) %T>% 
  saveRDS(., file.path(outpath, "aligned_UTR_sequences.mm.rn.bt.hs.RDS"))


#### MOUSE-HUMAN ORTHOLOGS ####
## get mouse-human 1to1 orthologs from ENSEMBL
# load Mart of mouse database from ensembl
# sometimes function useMart isn't able to connect to server and returns error, this chunck repeats useMart until there is no error
mart <- "error"
count <- 0
while(class(mart) == "character"){
  
  count <- count + 1
  print(count)
  
  # load ENSEMBL mart
  mart <- tryCatch(expr = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "http://jul2018.archive.ensembl.org"),
                   error = function(x) return("error"))
  
  # stop if count get too big
  if(count > 2){
    stop("Something's not right")
  }
  
}

# get info about genes, save
ensembl_info <- 
  getBM(attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_orthology_type"), mart = mart) %>% 
  as.tibble(.) %>% 
  dplyr::filter(hsapiens_homolog_orthology_type == "ortholog_one2one") %>% 
  dplyr::rename(gene_id = ensembl_gene_id) %T>% 
  readr::write_csv(., file.path(outpath, "ensembl.93.mouse_human.one2one_homologs.csv"))


##### 
# ## get sequences of mouse 3' UTR
# # get 3' UTR coordinates, split by gene and reduce ranges
# ensembl_3UTRs <-
#   ensembl_gtf %>%
#   gtfToGRanges(., filter = "three_prime_utr") %>% 
#   GenomicRanges::split(., .$gene_id) %>%
#   GenomicRanges::reduce(., ignore.strand = F) %>%
#   unlist(.)
# 
# # make names unique 
# names(ensembl_3UTRs) <- make.unique(names(ensembl_3UTRs))
# 
# # # save coordinates as .bed
# rtracklayer::export.bed(object = ensembl_3UTRs, con = file.path(outpath, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.3primeUTRs.bed"))
# 
# # get sequences of all 3' uTRs, save
# ensembl_3UTRs_seq <- 
#   BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, ensembl_3UTRs) %T>%
#   saveRDS(., "mm10.3pUTRs.DNAStringSet.RDS")


##### 
# ## clean aligned 3'UTRs
# # split mouse coordinates and sequences to 2 columns
# aligned_UTRs[, `:=`(mouse_coords = str_remove_all(mouse_coords_seq, "^mm10\\.|\\|.*"), 
#                     mouse_seq = str_remove(mouse_coords_seq, ".*\\|"), 
#                     mouse_coords_seq = NULL)]
# 
# # set column order
# setcolorder(aligned_UTRs, c("mouse_coords", "mouse_seq", "rn5_coords", "hg19_coords", "canFam3_coords", "rn5_seq", "hg19_seq", "canFam3_seq"))
# 
# # get mouse coordinates
# mouse_granges <- 
#   tibble(coords = aligned_UTRs$mouse_coords, 
#          mouse_coords = aligned_UTRs$mouse_coords) %>% 
#   tidyr::separate(., coords, c("seqnames", "start", "end"), sep = ":") %>% 
#   GRanges(.) 
# 
# # findOverlaps between aligned parts of 3' UTRs and whole 3'UTRs, get ENSEMBL gene IDs
# overlaps <- findOverlaps(mouse_granges, ensembl_3UTRs, type = "within", ignore.strand = T)
# mouse_overlaps <- mouse_granges[queryHits(overlaps)]
# mcols(mouse_overlaps)$gene_id <- names(ensembl_3UTRs[subjectHits(overlaps)])
# 
# # convert to data.table, join with aligned UTRs
# aligned_UTRs_mouse <- 
#   data.table(mouse_coords = mouse_overlaps$mouse_coords, 
#              gene_id = mouse_overlaps$gene_id) %>% 
#   .[aligned_UTRs, `:=`(mm10_seq = i.mouse_seq, 
#                        rn5_seq = i.rn5_seq, 
#                        hg19_seq = i.hg19_seq, 
#                        canFam3_seq = i.canFam3_seq), on = "mouse_coords"] %>% 
#   .[, `:=`(mouse_coords = str_c(mouse_coords, "|", gene_id),
#            gene_id = NULL)] %>% 
#   .[]
# 
# # create list of DNAStringSets
# aligned_3pUTRs_list <- list(mm10 = DNAStringSet(x = set_names(aligned_UTRs_mouse$mm10_seq, aligned_UTRs_mouse$mouse_coords)), 
#                             rn5 = DNAStringSet(x = set_names(aligned_UTRs_mouse$rn5_seq, aligned_UTRs_mouse$mouse_coords)), 
#                             hg19 = DNAStringSet(x = set_names(aligned_UTRs_mouse$hg19_seq, aligned_UTRs_mouse$mouse_coords)), 
#                             canFam3 = DNAStringSet(x = set_names(aligned_UTRs_mouse$canFam3_seq, aligned_UTRs_mouse$mouse_coords)))
# 
# # save RDS
# saveRDS(aligned_3pUTRs_list, file = file.path(outpath, "aligned_3pUTRs.mm10.rn5_hg19_canFam3.DNAStringSet.list.RDS"))


# # get ensembl-entrez ID mappings
# entrez_df <-
#   getBM(attributes = c("ensembl_gene_id", "entrezgene"), mart = mart) %>%
#   as.tibble(.) %>%
#   dplyr::rename(gene_id = ensembl_gene_id, entrez_id = entrezgene) %>%
#   dplyr::mutate(entrez_id = as.character(entrez_id))
# 
# 
# #####
# ## get genes and their miRNA target sites from targetScan
# # get list of genes
# genes <- ls(targetscan.Mm.egTARGETSFULL)
# 
# # get full info about miRNA targeting genes
# mirna_info_full <- mget(genes, targetscan.Mm.egTARGETSFULL)
# 
# # set slot names
# slot_names <- c("miRFamily", "UTRstart", "UTRend", "MSAstart", "MSAend", "Seedmatch", "PCT")
# 
# # reshape class targetscanTarget to tibble
# gene_mirna_targets <- purrr::map(names(mirna_info_full), function(gene_name){
# 
#   # go through slots, bind them to tibble
#   purrr::map(slot_names, function(slot){
# 
#     # get one slot
#     slot(mirna_info_full[[gene_name]], slot)
# 
#   }) %>%
#     bind_cols(.) %>%
#     set_colnames(., slot_names) %>%
#     dplyr::mutate(entrez_id = gene_name)
# 
# }) %>%
#   bind_rows(.) %>%
#   distinct(.) %>%
#   left_join(., entrez_df, by = "entrez_id") %T>%
#   readr::write_csv(., file.path(outpath, "targetScan", "targetScan.package.Mm.miRNA_gene_targets.csv"))

