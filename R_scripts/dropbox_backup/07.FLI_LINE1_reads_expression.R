### INFO: 
### DATE: Mon Aug 17 09:48:59 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.oocyte/map_clusters_to_Siomi/LINE1_FLI_expression")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(purrr)

library(GenomicRanges)
library(GenomicAlignments)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# library size path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.deduplicated.smallRNAseq/Data/Mapped/STAR_Siomi.multimappers"
library_size_path <- file.path(mapped_path, "4_library_size/library_sizes.txt")

# elements path
fli_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE"
fli_path <- file.path(fli_path, "LINE1.FLI_elements.bed")

# repeatMasker path
genome_path <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"
rmsk_clean_path <- file.path(genome_path, "rmsk.Siomi.20200701.clean.fa.out.gz")

# repeatModeler path
rmod_clean_path <- file.path(genome_path, "RepeatModeler/RepeatMasker")
rmod_clean_path <- file.path(rmod_clean_path, "rmsk.Siomi.20200728.RepeatModeler.clean.fa.out.gz")

# encode annotation path
encode_path <- file.path(genome_path, "ensembl.99.MesAur1.0.20200415.Siomi.UCSCseqnames.geneInfo.csv")

# mapped clusters path
bam_path <- file.path(inpath, "02_bam_LINE1_all_alignments")
bam_path <- list.files(bam_path, "*.bam$", full.names = T)

######################################################## READ DATA
# read library size table
library_size_tb <- readr::read_delim(library_size_path, delim = "\t", col_names = c("sample_id", "library_size"))

# read FLI elements coordinates
fli_gr <- rtracklayer::import(fli_path)

# read clean repeatMasker
rmsk_tb <- readr::read_delim(file = rmsk_clean_path, delim = "\t")

# read clean repeatModeler
rmod_tb <- readr::read_delim(file = rmod_clean_path, delim = "\t")
  
# read ENCODE annotation
gene_tb <- readr::read_csv(encode_path)

# get bam
bam_list <- purrr::map(bam_path, function(path){
  
  # reads with NH tag
  GenomicAlignments::readGAlignmentsList(path, use.names = T, param = ScanBamParam(tag = "NH")) %>% 
    unlist(.)
  
}) %>% 
  set_names(basename(bam_path) %>% str_remove(., "\\.bam"))

######################################################## MAIN CODE
# get LINE as GRanges
line_gr <- 
  rbind(rmsk_tb, rmod_tb) %>% 
  dplyr::filter(repClass %in% c("LINE"), 
                repFamily == "L1") %>% 
  # tidyr::unite(rmsk_id, seqnames, start, end, rmsk_id, repName, sep = ".", remove = F) %>% 
  GRanges(.)

# get ENCODE annotation as GRanges
gene_gr <- GRanges(gene_tb)

# get coverage
coverage_gr <-
  unname(bam_list) %>%
  do.call(c, .) %>%
  GRanges(.) %>%
  GenomicRanges::reduce(., ignore.strand = T, min.gapwidth = 10000)

# add ID
mcols(coverage_gr)$coverage_id <- str_c(seqnames(coverage_gr), start(coverage_gr), end(coverage_gr), sep = ".")

# summarize overlaps between reads and clusters
sum_overlaps <- purrr::map(names(bam_list), function(sample_id){
  
  # get one bam
  bam <- bam_list[[sample_id]]
  
  # find overlaps
  overlaps <- findOverlaps(coverage_gr, bam, ignore.strand = T)
  
  # for each cluster get count of reads, both full reads and fractions for multimappers
  count_tb <- 
    tibble(coverage_id = coverage_gr[queryHits(overlaps)]$coverage_id, 
           count = rep(1, length(subjectHits(overlaps)))) %>% 
    dplyr::group_by(coverage_id) %>% 
    dplyr::summarise(count = sum(count)) %>% 
    dplyr::mutate(sample_id = sample_id)
  
  # return 
  return(count_tb)
  
}) %>% 
  dplyr::bind_rows(.)

# get RPM
rpm_tb <-
  sum_overlaps %>% 
  dplyr::left_join(., library_size_tb, by = "sample_id") %>% 
  dplyr::mutate(library_size = (library_size / 1e6), 
                rpm = count / library_size) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "WT|KO")) %>% 
  dplyr::group_by(coverage_id, genotype) %>% 
  dplyr::summarise(rpm = mean(rpm)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::pivot_wider(id_cols = coverage_id, names_from = genotype, values_from = rpm, 
                     values_fill = 0, names_sep = ".", names_prefix = "RPM.") %>% 
  dplyr::arrange(-RPM.KO)

# get as GRanges
rpm_gr <- 
  rpm_tb %>% 
  dplyr::select(coverage_id) %>% 
  tidyr::separate(coverage_id, into = c("seqnames", "start", "end"), sep = "\\.", remove = F) %>% 
  GRanges(.)


### overlap with FLIs
# overlap
overlap <- findOverlaps(rpm_gr, fli_gr, ignore.strand = T)

# add overlap to the table
overlap_tb_1 <- tibble(coverage_id = rpm_gr[queryHits(overlap)]$coverage_id,
                     fli_hit_id = fli_gr[subjectHits(overlap)]$name)


### overlap with LINE1s
# overlap
overlap <- findOverlaps(rpm_gr, line_gr, ignore.strand = T)

# add overlap to the table
overlap_tb_2 <- tibble(coverage_id = rpm_gr[queryHits(overlap)]$coverage_id,
                       line1_id = line_gr[subjectHits(overlap)]$rmsk_id)


### overlap with ENCODE genes
# overlap
overlap <- findOverlaps(rpm_gr, gene_gr, ignore.strand = T)

# add overlap to the table
overlap_tb_3 <- tibble(coverage_id = rpm_gr[queryHits(overlap)]$coverage_id,
                       gene_id = gene_gr[subjectHits(overlap)]$gene_id, 
                       gene_name = gene_gr[subjectHits(overlap)]$gene_name)


### join expression table with annotation
# join with expression table
rpm_tb_annotated <- 
  rpm_tb %>% 
  dplyr::left_join(., overlap_tb_1, by = "coverage_id") %>% 
  dplyr::left_join(., overlap_tb_2, by = "coverage_id") %>% 
  dplyr::left_join(., overlap_tb_3, by = "coverage_id") %>% 
  dplyr::mutate(coordinates = str_replace_all(coverage_id, "\\.", " ")) %>% 
  dplyr::select(coordinates, RPM.WT, RPM.KO, fli_hit_id, line1_hit_id = line1_id, gene_id, gene_name) %>% 
  dplyr::mutate(fli_hit_id = replace(fli_hit_id, is.na(fli_hit_id), "no_hit"), 
                line1_hit_id = replace(line1_hit_id, is.na(line1_hit_id), "no_hit")) %>% 
  dplyr::group_by(coordinates) %>% 
  dplyr::slice_sample(n = 1) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::arrange(-RPM.KO)

# write
readr::write_csv(rpm_tb_annotated, file.path(outpath, "LINE1_FLI_reads.RPM.oocyte_deduplexed.csv"))

# write filtered
rpm_tb_annotated %>% 
  dplyr::filter(fli_hit_id == "no_hit", line1_hit_id == "no_hit") %T>% 
  readr::write_csv(., file.path(outpath, "LINE1_FLI_reads.RPM.oocyte_deduplexed.no_LINE1_hits.csv"))

