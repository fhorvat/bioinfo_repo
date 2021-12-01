### INFO: 
### DATE: Mon Oct 29 16:48:51 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids")

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

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "gtfToGRanges.R"))

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# miRBase gtf path
mirbase_gtf_path <- list.files(genome_path, pattern = "miRBase.*gff3", full.names = T)

# ensembl .gtf path
ensembl_path <- list.files(path = genome_path, pattern = "ensembl.93.*[0-9]{6}.UCSC.*gtf.gz", full.names = T)

# gene info path
gene_info_path <- list.files(path = genome_path, pattern = "ensembl.93.*[0-9]{6}.UCSC.*geneInfo.csv", full.names = T)

# miRBase mature miRNAs sequences path
mirbase_seq_path <- "/common/DB/genome_reference/miRBase/miRBase.22.20181029.mature.fa.gz"

# grid path
grid_path <- file.path(inpath, "grid.Fugaku.Encode.binned_fpkm.relative_order.100.csv")

######################################################## READ DATA
# read miRBase gff
mirna_gr <- rtracklayer::import.gff(con = mirbase_gtf_path)

# read ensembl .gtf
ensembl_gtf <- read_delim(file = ensembl_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# read gene info
genes_info <- readr::read_csv(gene_info_path)

# read miRBase sequences
mirbase_seq <- Biostrings::readRNAStringSet(filepath = mirbase_seq_path)

# read grid
grid_df <- readr::read_csv(grid_path)

######################################################## MAIN CODE
#####
## prepare data
# clean grid
grid_df %<>% 
  dplyr::select(gene_id, starts_with("bin"))

# prepare full grid
full_grid_df <- 
  as.tibble(expand.grid(1:100, 1:100)) %>% 
  set_colnames(c("bin_fpkm", "bin_fpkm_relative")) %>% 
  dplyr::mutate(bin_id = str_c(bin_fpkm, bin_fpkm_relative, sep = "."))

# prepare distance data.table
dist_dt <-
  full_grid_df %>% 
  column_to_rownames("bin_id") %>%
  dist(., method = "euclidean") %>% 
  as.matrix(.) %>%
  as.data.table(., keep.rownames = "bin_id") %>% 
  melt(., id.vars = "bin_id", variable.name = "bin_id_dist", value.name = "dist")

##### 
## 3' UTR
# get 3' UTR coordinates, split by gene and reduce ranges
ensembl_3UTRs <-
  ensembl_gtf %>%
  gtfToGRanges(., filter = "three_prime_utr") %>% 
  GenomicRanges::split(., .$gene_id) %>%
  GenomicRanges::reduce(., ignore.strand = F) %>%
  unlist(.)

# make names unique 
names(ensembl_3UTRs) <- make.unique(names(ensembl_3UTRs))

# get sequences of all 3' LTRs
ensembl_3UTRs_seq <- BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, ensembl_3UTRs)


##### 
## miRNA
# prepare ranges of mature miRNA
mirna_gr <- mirna_gr[mcols(mirna_gr)$type == "miRNA"]
mcols(mirna_gr) <- mcols(mirna_gr)[, c("Name")]
names(mcols(mirna_gr)) <- "gene_id"

# get miRNA coordinates
mirna_coords <- mirna_gr[str_detect(mcols(mirna_gr)$gene_id, "mmu-let-7d-5p")][1]

# get let7a sequences
# mirna_seq <- BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, mirna_coords)
mirna_seq <- 
  mirbase_seq[str_detect(string = names(mirbase_seq), pattern = mcols(mirna_coords)$gene_id)] %>% 
  Biostrings::DNAStringSet(.) %>% 
  unlist(.)

# find 3' UTR matching nucleotides 2-8
match_2to8 <- 
  Biostrings::vmatchPattern(pattern = subseq(mirna_seq, start = 2, end = 8) %>% reverseComplement(.), 
                            subject = ensembl_3UTRs_seq, 
                            max.mismatch = 0) %>% 
  unlist(.) %>% 
  as.tibble(.) %>% 
  dplyr::mutate(gene_id = str_remove(names, "\\..*$")) %>% 
  dplyr::count(gene_id, sort = T) %>% 
  dplyr::inner_join(., grid_df, by = "gene_id") %>% 
  dplyr::group_by(bin_fpkm, bin_fpkm_relative) %>% 
  dplyr::summarise(n = sum(n)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(bin_id = str_c(bin_fpkm, bin_fpkm_relative, sep = ".")) %>% 
  as.data.table(., key = bin_id) %>% 
  .[order(bin_fpkm, bin_fpkm_relative), .(n, bin_id)] %>% 
  .[dist_dt, on = c("bin_id" = "bin_id_dist")] %>% 
  .[is.na(n), n := 0] %>% 
  .[dist != 0] %>% 
  .[, dist := (1 / (dist ^ 2 + 0))] %>% 
  .[, n := (n * dist)] %>% 
  .[, .(n = sum(n)), by = i.bin_id] %>% 
  .[, c("bin_fpkm", "bin_fpkm_relative") := tstrsplit(i.bin_id, ".", fixed = TRUE)] %>% 
  .[, `:=`(bin_fpkm = as.integer(bin_fpkm), 
           bin_fpkm_relative = as.integer(bin_fpkm_relative), 
           n = as.vector(scale(n)))] %>% 
  .[order(bin_fpkm, bin_fpkm_relative)]

# match_2to8 <- 
#   match_2to8[n > 2, n := 2] %>% 
#   .[n < 0, n := 0] %>% 
#   .[]

# plot hits
ggplot(match_2to8, aes(bin_fpkm, bin_fpkm_relative)) + 
  geom_tile(aes(fill = n), colour = "white") + 
  scale_fill_gradientn(colours = c("blue", "green", "red"), 
                       values = scales::rescale(c(-1, 0, 1))) +
  ggsave(filename = file.path(outpath, "test6.png"))

