### INFO: read smallRNA seq bam file, get counts over exons
### DATE: Wed May 16 02:54:16 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# mouse
animal <- "mouse"
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/mouse/connected_regions")
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# # cow
# animal <- "cow"
# setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/cow/connected_regions")
# genome_path <- "/common/DB/genome_reference/cow/bosTau8.UMD3.1.GCA_000003055.4"

# # pig
# animal <- "pig"
# setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_2018May/smallRNA_expression/pig/connected_regions")
# genome_path <- "/common/DB/genome_reference/pig/susScr11.Sscrofa11.1.GCA_000003025.6"

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(purrr)
library(data.table)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "bamToGRangesList.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main outpath
outpath <- getwd()

# set main inpath
inpath <- getwd()

# tables list
table_list <- list.files(inpath, pattern = "smallRNA.connected_regions.*csv", full.names = T)

# gene info path
gene_info_path <- list.files(genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.geneInfo.csv", full.names = T)

# transcript info path
transcript_info_path <- list.files(genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.transcriptInfo.csv", full.names = T)

# get path of retroGenes
retro_path <- list.files(genome_path, pattern = "RetroGenesV6.*txt.gz", full.names = T)

# library size path
libsize_list <-
  list.files(file.path(inpath, "../../data_sets", animal), pattern = ".*lib_size.hist.txt",
             recursive = T, full.names = T) %>%
  tibble(lib_path = .) %>%
  dplyr::mutate(sample = basename(lib_path) %>% str_remove_all(., ".SE.lib_size.hist.txt|.SE_cl.lib_size.hist.txt|.SE.*"),
                experiment = str_remove(lib_path, "\\/Links.*") %>% basename(.) %>% str_remove(., "(?<=[0-9])_.*"))

######################################################## READ DATA
# read info about genes
genes_info <- 
  readr::read_csv(gene_info_path) %>% 
  tidyr::unite(coordinates_gene, seqnames, start, end, strand, sep = " ") %>% 
  dplyr::mutate(gene_description = stringr::str_remove(gene_description, " \\[.*"))

# read info about transcripts
transcripts_info <- readr::read_csv(transcript_info_path)

# read library sizes to one data.frame
retro_df <- readr::read_delim(file = retro_path, delim = "\t")

# read library sizes to one data.frame
libsize_df <-
  purrr::pmap(libsize_list, function(lib_path, sample, experiment) {
    
    # read library data.frame from files defined in libsize_list data.frame
    lib_df <-
      data.table::fread(lib_path) %>%
      as.tibble(.) %>%
      magrittr::set_colnames(c("lib_size", "seq_length")) %>%
      dplyr::mutate(ID = str_c(experiment, sample, sep = "."))
    
  }) %>%
  dplyr::bind_rows(.) %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(lib_size_all = sum(lib_size),
                   lib_size_19_24 = sum(lib_size[seq_length >= 19 & seq_length <= 24]),
                   lib_size_21_23 = sum(lib_size[seq_length >= 21 & seq_length <= 23])) %>%
  dplyr::mutate_at(vars(starts_with("lib_size")), funs(round(. / 1E6, 5)))

######################################################## MAIN CODE
# make transcript ID's unique
retro_df <-  
  dplyr::bind_rows(retro_df %>% dplyr::filter(mgc != "noMgc"), 
                   retro_df %>% dplyr::filter(mgc == "noMgc") %>% dplyr::mutate(mgc = make.unique(mgc)))

# read all tables, join with all categories
count_df <-
  purrr::map(table_list, readr::read_csv) %>%
  dplyr::bind_rows(.) %>%
  dplyr::filter((experiment != "Tam_2008" | sample == "s_oocyte_19to24")) %>%
  tidyr::unite(col = "ID", experiment, sample, sep = ".")

### retroGenes ranges
# get ranges of pseudogenes
retro_pseudo_gr <- 
  GenomicRanges::GRanges(seqnames = retro_df$`#chrom`, 
                         ranges = IRanges(start = retro_df$chromStart, end = retro_df$chromEnd), 
                         strand = retro_df$strand, 
                         pseudogene_id = retro_df$name) %>% 
  unique(.)

# get ranges of original transcripts which have annotated pseudogenes
retro_transcripts_gr <- 
  GenomicRanges::GRanges(seqnames = retro_df$gChrom, 
                         ranges = IRanges(start = retro_df$gStart, end = retro_df$gEnd), 
                         strand = retro_df$gStrand, 
                         gene_id = retro_df$mgc) %>% 
  unique(.) %>% 
  GenomicRanges::split(., .$gene_id) %>% 
  GenomicRanges::reduce(.) %>% 
  unlist(.)

mcols(retro_transcripts_gr)$gene_id <- names(retro_transcripts_gr)
names(retro_transcripts_gr) <- NULL

### connected regions ranges - separate ranges
# ranges 1 in pair of connected regions
regions_gr1 <- 
  count_df %>% 
  tidyr::separate(col = coordinates_origin, into = c("seqnames", "start", "end_strand"), sep = " ") %>% 
  tidyr::separate(col = end_strand, into = c("end", "strand"), sep = "\\|") %>% 
  dplyr::mutate(start = as.numeric(start), end = as.numeric(end)) %>% 
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)
colnames(mcols(regions_gr1))[colnames(mcols(regions_gr1)) == "coordinates_target"] <- "coordinates_connected"

# ranges 2 in pair of connected regions
regions_gr2 <- 
  count_df %>% 
  tidyr::separate(col = coordinates_target, into = c("seqnames", "start", "end_strand"), sep = " ") %>% 
  tidyr::separate(col = end_strand, into = c("end", "strand"), sep = "\\|") %>% 
  dplyr::mutate(start = as.numeric(start), end = as.numeric(end)) %>% 
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)
colnames(mcols(regions_gr2))[colnames(mcols(regions_gr2)) == "coordinates_origin"] <- "coordinates_connected"


### get overlaps with transcripts which have known pseudogenes
# overlap
transcripts_overlaps1 <- GenomicRanges::findOverlaps(regions_gr1, retro_transcripts_gr)
transcripts_overlaps2 <- GenomicRanges::findOverlaps(regions_gr2, retro_transcripts_gr)

# extract transcripts which share reads
connected_transcripts <- 
  dplyr::bind_rows(regions_gr1[queryHits(transcripts_overlaps1)] %>% as.data.frame(.), 
                   regions_gr2[queryHits(transcripts_overlaps2)] %>% as.data.frame(.)) %>% 
  as.tibble(.) %>% 
  dplyr::mutate(gene_id = retro_transcripts_gr[c(subjectHits(transcripts_overlaps1), subjectHits(transcripts_overlaps2))]$gene_id) %>% 
  unique(.)


### check if transcripts share reads with pseudogenes
# get connected regions 
pseudo_gr <-
  connected_transcripts %>% 
  dplyr::select(coordinates_connected, reg_con) %>% 
  tidyr::separate(col = coordinates_connected, into = c("seqnames", "start", "end_strand"), sep = " ") %>% 
  tidyr::separate(col = end_strand, into = c("end", "strand"), sep = "\\|") %>% 
  dplyr::mutate(start = as.numeric(start), end = as.numeric(end)) %>% 
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

# overlap with pseudogene annotation
pseudo_overlaps <- GenomicRanges::findOverlaps(pseudo_gr, retro_pseudo_gr)

# get all IDs of pseudogenes for each connected region  
pseudo_connected <- 
  dplyr::bind_cols(reg_con = pseudo_gr[queryHits(pseudo_overlaps)] %$% reg_con, 
                   retro_pseudo_gr[subjectHits(pseudo_overlaps)] %>% as.data.frame(.)) %>% 
  tidyr::unite(col = "coordinates_pseudogene", seqnames, start, end, strand, sep = " ") %>% 
  dplyr::select(-width)

### get transcript - pseudogene connected pairs
# create data.frame 
connected_all <- 
  
  # select columns from connected transcripts table
  connected_transcripts %>% 
  dplyr::select(reg_con, read_connection_count, ID, gene_id) %>% 
  
  # add pseudogene coordinates
  dplyr::right_join(., pseudo_connected, by = "reg_con") %>% 
  
  # summarize counts of reads in gene-pseudogene pairs
  dplyr::group_by(gene_id, pseudogene_id, ID) %>% 
  dplyr::summarise(read_connection_count = sum(read_connection_count)) %>% 
  dplyr::ungroup(.) %>% 
  
  # join with library sizes, normalize
  dplyr::left_join(., libsize_df %>% dplyr::select(ID, lib_size_21_23), by = "ID") %>% 
  dplyr::mutate(read_connection_fpm = round((read_connection_count / lib_size_21_23), 3)) %>% 
  dplyr::select(-lib_size_21_23) %>% 
  
  # gather, unite, spread
  tidyr::gather(variable, value, -(c(gene_id:ID))) %>% 
  dplyr::mutate(variable = str_remove(variable, pattern = "read_connection_")) %>% 
  tidyr::unite(ID, ID, variable, sep = ".") %>% 
  tidyr::spread(key = ID, value = value) %>% 
  
  # add pseudogene coordindates
  dplyr::left_join(., retro_pseudo_gr %>% as.data.frame(.), by = "pseudogene_id") %>% 
  tidyr::unite(col = coordinates_pseudo, seqnames, start, end, strand, sep = " ") %>% 
  
  # replace transcript coordinates with gene coordinates and add info about genes
  dplyr::mutate(gene_id = str_remove(gene_id, "\\..*")) %>% 
  dplyr::rename(transcript_id = gene_id) %>% 
  dplyr::left_join(transcripts_info, by = "transcript_id") %>% 
  dplyr::left_join(genes_info, by = "gene_id") %>%

  # replace NA with 0, calculate average values
  dplyr::mutate_at(vars(matches("count|fpm")), funs(replace(., is.na(.), 0))) %>% 
  dplyr::mutate(avg.count = rowMeans(select_at(., vars(contains("count")))), 
                avg.fpm = rowMeans(select_at(., vars(contains("fpm"))))) %>% 
  
  # tidy and arange
  dplyr::select(coordinates_gene, gene_id, gene_name, gene_description, gene_biotype, coordinates_pseudo, pseudogene_id, avg.count, contains("count"), avg.fpm, contains("fpm")) %>% 
  dplyr::arrange(desc(avg.fpm))


  
  
