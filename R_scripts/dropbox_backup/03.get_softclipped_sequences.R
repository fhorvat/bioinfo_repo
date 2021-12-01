### INFO: 
### DATE: Tue Aug 03 14:21:57 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/functional_oocyte_miRNA/maternal_miRNA_expression/nontemplated_additions/bam_subset")

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
library(rtracklayer)
library(GenomicAlignments)
library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# list of miRNA info files
mirna_info_list <- list.files(file.path(inpath, ".."), ".*\\.geneInfo.csv$", full.names = T)

# list of bam files
bam_list <- 
  list.files(inpath, ".*\\.bam$", full.names = T) %>% 
  .[!str_detect(., "s_oocyte_large.21to23nt")]

######################################################## READ DATA
# read miRNA info
mirna_info <- purrr::map(mirna_info_list, function(path){
  
  readr::read_csv(path) %>% 
    dplyr::mutate(animal = str_extract(path, "mouse.mm10.GarciaLopez_2015_RNA_GSE59254|cow.bosTau9|pig.susScr11"))
    
})

# read bam files
mirna_bam_list <- purrr::map(bam_list, function(path){
  
  # read bam
  bam_gr <- GenomicAlignments::readGAlignmentsList(file = path, 
                                                   param = ScanBamParam(what = c("qname", "seq"), 
                                                                        tag = c("NH"), 
                                                                        flag = scanBamFlag(isSecondaryAlignment = NA))) 
  
  # return
  return(bam_gr)
  
}) %>% 
  set_names(., str_extract(bam_list, "mouse.mm10.GarciaLopez_2015_RNA_GSE59254|cow.bosTau9|pig.susScr11"))

######################################################## MAIN CODE
# set dataset names
dataset_names <- c("mouse.mm10.GarciaLopez_2015_RNA_GSE59254", "cow.bosTau9", "pig.susScr11")

# create miRNA info table
mirna_info_tb <- 
  mirna_info %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::filter(!(gene_id %in% c("MIMAT0013865_1", "MIMAT0002152_1", "MIMAT0007754_1", "MIMAT0046670", "MIMAT0016938", 
                                 "MIMAT0003516_1", "MIMAT0009383_1", 
                                 "MIMAT0000525"))) %>%
  dplyr::group_by(animal) %>% 
  dplyr::top_n(5, FPM_oocyte) %>%
  dplyr::ungroup(.)
  
# prepare miRNA coordinates and names in table
mirna_coordinates_tb <- 
  mirna_info_tb %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::mutate(mirna_start = ifelse(strand == "+", start, end)) %>% 
  dplyr::select(gene_id, mirna_start) 

# get GRanges
mirna_bed_list <- 
  mirna_info_tb %>% 
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.) %>% 
  split(., mcols(.)$animal)

# assign reads to miRNAs
mirna_reads <- purrr::map(dataset_names, function(animal){
  
  # get bam
  animal_bam <- 
    mirna_bam_list[[animal]] %>% 
    unlist(.)
  
  # get bed
  animal_bed <- mirna_bed_list[[animal]]
  
  # overlap
  overlaps <- findOverlaps(animal_bam, animal_bed, ignore.strand = F)
  
  # assing miRNA to each read
  animal_mirna <- animal_bam[queryHits(overlaps)]
  mcols(animal_mirna)$gene_id <- mcols(animal_bed[subjectHits(overlaps)])$gene_id
  mcols(animal_mirna)$animal <- animal
  
  # return
  return(animal_mirna)
  
}) %>% 
  set_names(dataset_names)


# take only reads which have the same start as annotated miRNA 
# for + strand miRNAs this is "start" in GRanges, for - strand miRNAs it's the "end" in GRanges
# get reverse complement sequences of miRNAs on minus strand
# get soft clipped cigars 
# for + strand miRNAs it's "S" in the end of cigar string, for - strand miRNAs it's the "S" in beginning of cigar string
# remove reads with 5' end softclipping
# extract read sequence and softclip sequence
mirna_sequences <- 
  purrr::map(mirna_reads, as_tibble) %>%
  dplyr::bind_rows(.) %>% 
  dplyr::select(seqnames, start, end, strand, seq, gene_id, NH, cigar, animal) %>% 
  dplyr::mutate(read_start = ifelse(strand == "+", start, end)) %>% 
  dplyr::left_join(., mirna_coordinates_tb, by = "gene_id") %>% 
  dplyr::filter(read_start == mirna_start) %>% 
  dplyr::filter(!str_detect(seq, "N")) %>% 
  dplyr::mutate(seq_revcomp = ifelse(strand == "-", DNAStringSet(seq) %>% reverseComplement(.) %>% as.character(.), seq), 
                soft_clip_start = str_extract(cigar, "^[0-9]+S(?=[0-9]{2}M)"), 
                soft_clip_end = str_extract(cigar, "(?<=[0-9]{2}M)[0-9]+S$")) %>% 
  dplyr::filter((strand == "+" & is.na(soft_clip_start)) | (strand == "-" & is.na(soft_clip_end))) %>%
  dplyr::mutate(cigar_softclip = ifelse(is.na(soft_clip_start), soft_clip_end, soft_clip_start), 
                cigar_match = str_remove(cigar, cigar_softclip)) %>% 
  dplyr::select(seqnames:strand, cigar, seq = seq_revcomp, cigar_match, cigar_softclip, NH, gene_id, animal) %>% 
  dplyr::mutate(read_seq = str_extract(seq, str_c("^.{", str_remove(cigar_match, "M"), "}")), 
                softclip_seq = str_remove(seq, read_seq))

# get miRNA expression
mirna_expression <-
  mirna_sequences %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::summarise(mirna_count = sum(1 / NH))

# # remove reads from Garcia-Lopez which have softclipped "ATC" sequence (adapter contamination)
# mirna_sequences %<>% 
#   dplyr::filter(animal != "mouse.mm10.GarciaLopez_2015_RNA_GSE59254" | !(str_detect(softclip_seq, "^ATC")) | is.na(softclip_seq))


### plot frequency of each nucleotide post 3' end, separate for each animal
# create polyN sequences in a table
nontemplated_additions <- purrr::map(c("A", "U", "G", "C"), function(nucleotide){
  
  purrr::map(1:4, function(n) rep(nucleotide, n) %>% unlist(.) %>% str_c(., collapse = ""))
  
}) %>% 
  unlist(.) %>% 
  tibble(softclip_seq = .) %>% 
  dplyr::mutate(nucleotide = str_extract(softclip_seq, "A|U|C|G"), 
                nucleotide = factor(nucleotide, levels = c("A", "U", "G", "C"))) 
  
# get counts of sofclipped reads per animal and gene
softclipped_counts <- 
  mirna_sequences %>% 
  dplyr::filter(softclip_seq %in% nontemplated_additions$softclip_seq) %>% 
  dplyr::group_by(animal, gene_id) %>% 
  dplyr::summarize(nta_count = n()) %>% 
  dplyr::left_join(., mirna_expression, by = "gene_id") %>% 
  dplyr::left_join(., mirna_info_tb, by = c("animal", "gene_id")) %>% 
  dplyr::arrange(animal, desc(FPM_oocyte))
  
  
# separate table to animals
mirna_animals <- split(mirna_sequences, mirna_sequences$animal)

# count sequences post 3' end for each animal
mirna_animals_counts <- purrr::map(names(mirna_animals), function(animal){
  
  # get counts of different sequences
  mirna_counts_tb <- 
    mirna_animals[[animal]] %>%
    dplyr::filter(!is.na(softclip_seq)) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::count(softclip_seq) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(gene_id, desc(n))

  # get only polyN
  mirna_poly <- 
    mirna_counts_tb %>% 
    dplyr::mutate(softclip_seq = str_replace(softclip_seq, "T", "U")) %>% 
    dplyr::inner_join(., nontemplated_additions, by = "softclip_seq") %>% 
    dplyr::left_join(., mirna_info_tb, by = "gene_id") %>%
    dplyr::left_join(., mirna_expression, by = "gene_id") %>% 
    dplyr::mutate(percentage = 100 * (n / mirna_count)) %>% 
    dplyr::mutate(softclip_seq = factor(softclip_seq, levels = rev(nontemplated_additions$softclip_seq))) %>% 
    dplyr::arrange(-FPM_oocyte) %>% 
    dplyr::mutate(gene_name = factor(gene_name, levels = unique(gene_name))) %>% 
    dplyr::select(gene_name, nucleotide, softclip_seq, percentage, nta_count = n, mirna_count)

  ### save as wide table
  # prepare table
  mirna_wide <- 
    mirna_poly %>% 
    dplyr::select(gene_name, nta_seq = softclip_seq, nta_count, mirna_count, percentage)
  
  # get the count table
  mirna_wide_count <- 
    mirna_wide %>% 
    dplyr::mutate(nta_seq = forcats::fct_rev(nta_seq)) %>% 
    complete(nta_seq = nta_seq) %>% 
    tidyr::pivot_wider(id_cols = c(gene_name, mirna_count), names_from = nta_seq, values_from = nta_count) %>% 
    dplyr::filter(!is.na(gene_name)) %>% 
    replace(is.na(.), 0) %T>%
    readr::write_csv(., file = file.path(outpath, str_c("miRNA_3p_NTA.barplot", animal, "csv", sep = ".")))
  
  # plot
  nta_barplot <- 
    ggplot() + 
    geom_bar(data = mirna_poly, 
             mapping = aes(x = nucleotide, y = percentage, fill = softclip_seq), 
             width = 1, stat = "identity", position = "stack", color = "black") + 
    facet_wrap(~ gene_name, nrow = 1, strip.position = "bottom") +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(0, 60)) + 
    scale_fill_manual(values = c("A" = "#993333", "AA" = "#B76D6D", "AAA" = "#D4A8A8", "AAAA" = "#F3E3E3", 
                                 "U" = "#4166AF", "UU" = "#738EC5", "UUU" = "#A6B7DC", "UUUU" = "#D9E0F3", 
                                 "G" = "#7B6BB2", "GG" = "#9D91C6", "GGG" = "#BFB7DA", "GGGG" = "#E2DEEF",
                                 "C" = "#6F8940", "CC" = "#94A771", "CCC" = "#BAC5A3", "CCCC" = "#E0E4D5"), 
                      breaks = c("A", "AA", "AAA", "AAAA",
                                 "U", "UU", "UUU", "UUUU",
                                 "G", "GG", "GGG", "GGGG",
                                 "C", "CC", "CCC", "CCCC"),
                      labels = c("A", "2A", "3A", "4A",
                                 "U", "2U", "3U", "4U",
                                 "G", "2G", "3G", "4G",
                                 "C", "2C", "3C", "4C"),
                      drop = F) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(), 
          legend.title = element_blank()) +
    # theme(panel.border = element_blank()) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom", 
          strip.background = element_blank())

  # save
  ggsave(filename = file.path(outpath, str_c("miRNA_3p_NTA.barplot", animal, "png", sep = ".")),
         plot = nta_barplot,
         width = 12,
         height = 10,
         limitsize = FALSE)
  
  # return
  return(animal)
  
})


