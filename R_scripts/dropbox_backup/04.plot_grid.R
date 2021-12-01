### INFO: 
### DATE: Sun Nov 04 19:08:17 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/Su_2004_ProcNatlAcadSciUSA_GSE1133")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

library(ggthemes)

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
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# miRBase gtf path
mirbase_gtf_path <- list.files(genome_path, pattern = "miRBase.*gff3", full.names = T)

# ensembl .gtf path
ensembl_path <- list.files(path = genome_path, pattern = "ensembl.93.*[0-9]{6}.UCSC.*gtf.gz", full.names = T)

# miRBase mature miRNAs sequences path
mirbase_seq_path <- "/common/DB/genome_reference/miRBase/miRBase.22.20181029.mature.fa.gz"

# grid path
grid_path <- file.path(inpath, "grid.Su.mas5.oocyte.61.bins.csv")

# multiple alignment of 3' UTRs between mouse, rat, dog and human 
aligned_3pUTRs_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/multiz60way.3primeUTRs/aligned_3pUTRs.mm10.rn5_hg19_canFam3.RDS" 

######################################################## READ DATA
# read miRBase gff
mirna_gr <- rtracklayer::import.gff(con = mirbase_gtf_path)

# read ensembl .gtf
ensembl_gtf <- read_delim(file = ensembl_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# read miRBase sequences
mirbase_seq <- Biostrings::readRNAStringSet(filepath = mirbase_seq_path)

# read and clean grid
grid_df <- 
  readr::read_csv(grid_path) %>% 
  dplyr::select(gene_id, starts_with("bin"))

# read multiple alignment of 3' UTRs between mouse, rat, dog and human  
aligned_3pUTRs <- readRDS(aligned_3pUTRs_path)

######################################################## MAIN CODE
#####
## prepare data
# prepare full grid
full_grid_df <- 
  as.tibble(expand.grid(1:61, 1:61)) %>% 
  set_colnames(c("bin_absolute", "bin_relative")) %>% 
  dplyr::mutate(bin_id = str_c(bin_absolute, bin_relative, sep = "."))

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
# rtracklayer::export.bed(object = ensembl_3UTRs, con = file.path(outpath, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.3primeUTRs.bed"))

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
  mirbase_seq[str_detect(string = names(mirbase_seq), pattern = str_c(mcols(mirna_coords)$gene_id, collapse = "|"))] %>% 
  Biostrings::DNAStringSet(.) %>% 
  unlist(.)


#####
## aligned 3' UTRs
aligned_3pUTRs_mm10 <- aligned_3pUTRs[str_detect(names(aligned_3pUTRs), "^mm10")]
aligned_3pUTRs_rn5 <- aligned_3pUTRs[str_detect(names(aligned_3pUTRs), "^rn5")]
aligned_3pUTRs_canFam3 <- aligned_3pUTRs[str_detect(names(aligned_3pUTRs), "^canFam3")]
aligned_3pUTRs_hg19 <- aligned_3pUTRs[str_detect(names(aligned_3pUTRs), "^hg19")]

##### 
## create and plot grid
# set smoothing factor k
k <- 10

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
  dplyr::group_by(bin_absolute, bin_relative) %>% 
  dplyr::summarise(n = sum(n)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(bin_id = str_c(bin_absolute, bin_relative, sep = ".")) %>% 
  as.data.table(., key = bin_id) %>% 
  .[order(bin_absolute, bin_relative), .(n, bin_id)] %>% 
  .[dist_dt, on = c("bin_id" = "bin_id_dist")] %>% 
  .[is.na(n), n := 0] %>% 
  .[dist != 0]

# smooth 
match_2to8_smooth <- 
  copy(match_2to8) %>% 
  .[, dist := (1 / (dist ^ 2 + k))] %>% 
  .[, n := (n * dist)] %>% 
  .[, .(n = sum(n)), by = i.bin_id] %>% 
  .[, n := as.vector(scale(n))] %>% 
  .[, c("bin_absolute", "bin_relative") := lapply(tstrsplit(i.bin_id, ".", fixed = TRUE), as.integer)] %>% 
  .[order(bin_absolute, bin_relative)]

# plot hits
ggplot(match_2to8_smooth, aes(bin_absolute, bin_relative)) + 
  geom_tile(aes(fill = n)) + 
  scale_fill_gradientn(colours = c("blue", "green", "red"),
                       values = scales::rescale(c(-1, 0, 1))) +
  theme_tufte() +
  ggsave(filename = file.path(outpath, str_c("grid.Su.mas5.oocyte.61.bins.mmu-let-7-5p.", k, ".smooth.png")))


