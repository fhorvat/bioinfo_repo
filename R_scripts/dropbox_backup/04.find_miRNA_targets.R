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
library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# grid path
grid_path <- file.path(inpath, "grid.Su.mas5.skeletalmuscle.61.bins.csv")

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# miRBase gtf path
mirbase_gtf_path <- list.files(genome_path, pattern = "miRBase.*gff3", full.names = T)

# miRBase mature miRNAs sequences path
mirbase_seq_path <- "/common/DB/genome_reference/miRBase/miRBase.22.20181029.mature.fa.gz"

# accessory data path
accessory_data_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/accessory_data"

# mouse 3' UTR sequences path
mouse_3UTRs_path <- file.path(accessory_data_path, "mm10.3pUTRs.DNAStringSet.RDS")

# multiple alignment of 3' UTRs between mouse, rat, dog and human 
aligned_3pUTRs_path <- file.path(accessory_data_path, "aligned_3pUTRs.mm10.rn5_hg19_canFam3.DNAStringSet.list.RDS")

# distance grid
distance_grid_path <- file.path(accessory_data_path, "distance_grid.Euclidean.61x61.dt.RDS")

######################################################## READ DATA
# read and clean grid
grid_df <- 
  readr::read_csv(grid_path) %>% 
  dplyr::select(gene_id, starts_with("bin"))

# read miRBase gff
mirna_gr <- rtracklayer::import.gff(con = mirbase_gtf_path)

# read miRBase sequences
mirbase_seq <- Biostrings::readRNAStringSet(filepath = mirbase_seq_path)

# read mouse 3' UTR sequences
mouse_3UTRs_seq <- readRDS(mouse_3UTRs_path)

# read multiple alignment of 3' UTRs between mouse, rat, dog and human  
aligned_3pUTRs <- readRDS(aligned_3pUTRs_path)

# read distance data.table
dist_dt <- readRDS(distance_grid_path)

######################################################## MAIN CODE
##### 
## prepare data
# prepare ranges of mature miRNA
mirna_gr <- mirna_gr[mcols(mirna_gr)$type == "miRNA"]
mcols(mirna_gr) <- mcols(mirna_gr)[, c("Name")]
names(mcols(mirna_gr)) <- "gene_id"

# remove gaps from aligned 3'UTRs
aligned_3pUTRs <- purrr::map(aligned_3pUTRs, DECIPHER::RemoveGaps)

##### 
## get miRNA targets in 3' UTRs
# get miRNA coordinates
mirna_coords <- mirna_gr[str_detect(mcols(mirna_gr)$gene_id, "mmu-miR-133[:alpha:]{1}-3p")]

# get chosen miRNA sequence
mirna_seq <- 
  mirbase_seq[str_detect(string = names(mirbase_seq), pattern = str_c(mcols(mirna_coords)$gene_id, collapse = "|"))] %>% 
  Biostrings::DNAStringSet(.) %>% 
  .[1] %>% 
  unlist(.)

# get seed sequence
seed_2to8 <- 
  subseq(mirna_seq, start = 2, end = 8) %>% 
  reverseComplement(.)

# find seed in aligned 3' UTRs = conserved miRNA target sites
match_2to8_conserved <- purrr::map(aligned_3pUTRs, function(aligned_3pUTR_animal){
  
  # match pattern in 3'UTR
  Biostrings::vmatchPattern(pattern = seed_2to8, 
                            subject = aligned_3pUTR_animal, 
                            max.mismatch = 0) %>% 
    unlist(.) %>% 
    as.data.table(.) %>% 
    setkey(., names)
  
}) %>% 
  purrr::reduce(., merge, all = T) %>% 
  .[complete.cases(.)] %>% 
  .[, "names"]

# find seed in mouse 3'UTRs = all miRNA target sites
match_2to8_all <- 
  Biostrings::vmatchPattern(pattern = seed_2to8, 
                            subject = mouse_3UTRs_seq, 
                            max.mismatch = 0) %>% 
  unlist(.) 
  

##### 
## create and plot grid
# find 3' UTR matching nucleotides 2-8
mirna_targets <- 
  match_2to8_all %>% 
  as.tibble(.) %>% 
  dplyr::mutate(gene_id = str_remove_all(names, "^chr.*\\||\\..*$")) %>% 
  dplyr::count(gene_id, sort = T) %>% 
  dplyr::inner_join(., grid_df, by = "gene_id") %>% 
  dplyr::group_by(bin_absolute, bin_relative) %>% 
  dplyr::summarise(n = sum(n)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(bin_id = str_c(bin_absolute, bin_relative, sep = "."))

# plot just scatterplot 
ggplot(mirna_targets, aes(bin_absolute, bin_relative)) +
  geom_point(aes(size = n)) +
  theme_bw() +
  ggsave(filename = file.path(outpath, "grid.scatter.Su.mas5.skeletalmuscle.61.bins.mmu-miR-133-3p.all.png"))

# plot heatmap
mirna_targets_dt <-
  mirna_targets %>% 
  as.data.table(., key = bin_id) %>% 
  .[order(bin_absolute, bin_relative), .(n, bin_id)] %>% 
  .[dist_dt, on = c("bin_id" = "bin_id_dist")] %>% 
  .[is.na(n), n := 0] %>% 
  .[]

# set smoothing factor k
k <- 30

# smooth 
mirna_targets_smooth <- 
  copy(mirna_targets_dt) %>% 
  .[, dist := (1 / (dist ^ 2 + k))] %>% 
  .[, n := (n * dist)] %>% 
  .[, .(n = sum(n)), by = i.bin_id] %>% 
  .[, n := as.vector(scale(n))] %>% 
  .[, c("bin_absolute", "bin_relative") := lapply(tstrsplit(i.bin_id, ".", fixed = TRUE), as.integer)] %>% 
  .[order(bin_absolute, bin_relative)]

# plot hits
ggplot(mirna_targets_smooth, aes(bin_absolute, bin_relative)) + 
  geom_tile(aes(fill = n)) + 
  scale_fill_gradientn(colours = c("blue", "green", "red"),
                       values = scales::rescale(c(-1, 0, 1))) +
  theme_tufte() +
  ggsave(filename = file.path(outpath, str_c("grid.heatmap.Su.mas5.skeletalmuscle.61.bins.mmu-miR-133-3p.", k, ".smooth.all.png")))


