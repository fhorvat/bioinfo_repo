### INFO: read coordinates, get sequences from ENSEMBL biomart
### DATE: 06. 10. 2017.
### AUTHOR: Vedran Franke, Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Shubhangini/Ago2_repeats_evolution/Sequences")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

library(biomaRt)
library(Biostrings)

######################################################## PATH VARIABLES
outpath <- getwd()

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## READ DATA
# set animals from which you want sequences 
animals_list <- c("mouse" = "mmusculus",
                  "rat" = "rnorvegicus",
                  "golden hamster" = "mauratus",
                  "chinese hamster" = "ccrigri", 
                  "kangaroo rat" = "dordii",
                  "naked mole rat" = "hmale",	
                  "guinea pig" = "cporcellus",
                  "squirrle" = "itridecemlineatus",
                  "rabbit" = "ocuniculus",
                  "cow" = "btaurus",
                  "human" = "hsapiens"
                  # "pig" = "sscrofa"
                  # "elephant" = "lafricana",
                  # "bat" = "mlucifugus",
                  # "platypus" = "oanatinus",
                  # "dog" = "cfamiliaris"
)

# read Ago2 coordinates of sequences 
coords <- readr::read_csv(file = file.path(outpath, "Ago2_repeats_regions_coords.csv"))

######################################################## MAIN CODE
# set geneID
gene_symbol <- "Ago2"

# clean loaded coordinates of sequences
coords_clean <- 
  coords %>% 
  tidyr::separate(col = "coordinates", into = c("chromosome", "coordinates"), sep = ":") %>% 
  tidyr::separate(col = "coordinates", sep = "-", into = c("start", "end")) %>% 
  dplyr::left_join(., tibble(animal = animals_list, animal_name = names(animals_list)), by = "animal") %>% 
  dplyr::mutate(dataset = str_c(animal, "_gene_ensembl", sep = ""))

# loops over the organisms and gets coordinates/sequences
seq_data <- 
  lapply(unique(coords_clean$dataset), function(ensembl_dataset){
    
    ensembl_dataset <- unique(coords_clean$dataset)[1]
    
    # filter coordinates data.frame
    coords_filt <- 
      coords_clean %>% 
      dplyr::filter(dataset == ensembl_dataset)
    
    # load Mart
    mart <- "error"
    count <- 0
    while(class(mart) == "character"){
      count <- count + 1
      print(str_c(ensembl_dataset, " ", count))
      mart <- tryCatch(expr = useMart(biomart =  "ensembl", dataset  = ensembl_dataset),
                       error = function(x) return("error"))
    }
    
    
    
    # loop over all features from same mart and get sequences
    mart_seq <- 
      lapply(unique(coords_filt$intron), function(feature){
        
        feature <- unique(coords_filt$intron)[2]
        
        # go over features one by one and filter
        feature_filt <- 
          coords_filt %>% 
          dplyr::filter(intron == feature)
        
        # get sequences
        feature_seq <- 
          biomaRt::getSequence(chromosome = feature_filt$chromosome, 
                               start = feature_filt$start, 
                               end = feature_filt$end,
                               seqType = "gene_flank", 
                               type = "ensembl_gene_id",
                               mart = mart) %>% 
          dplyr::mutate(feature = feature)
        
      }) %>% 
      dplyr::bind_rows(.)
      
    
    # 
    # # get sequence from mart
    # sequence <- getSequence(id = ids[str_replace(string = ensembl_dataset, pattern = "_gene_ensembl", replacement = "")], 
    #                         type = "ensembl_gene_id",
    #                         seqType = "gene_exon_intron", 
    #                         mart = mart)
    # 
    # # get ensembl_transcript_id of transcript which spans longest gene coordinates
    # transcript_id <- 
    #   getBM(attributes = c("chromosome_name", "transcript_start", "transcript_end", "strand", "ensembl_transcript_id"),
    #         filters = c("ensembl_gene_id"),
    #         values = ids[str_replace(string = ensembl_dataset, pattern = "_gene_ensembl", replacement = "")],
    #         mart = mart) %>% 
    #   dplyr::mutate(transcript_length = abs(transcript_start - transcript_end), 
    #                 dataset = str_replace(string = ensembl_dataset, pattern = "_gene_ensembl", replacement = "")) %>% 
    #   dplyr::top_n(n = 1, transcript_length) %$% 
    #   ensembl_transcript_id
    # 
    # # get coordinates and whole sequences of genes
    # coding_seq <-
    #   getBM(attributes = c("start_position", "end_position", "cdna"),
    #         filters = c("ensembl_transcript_id"),
    #         values = transcript_id,
    #         mart = mart) %>%
    #   dplyr::mutate(dataset = str_replace(string = ensembl_dataset, pattern = "_gene_ensembl", replacement = ""))
    # 
    # # get cooridnates of exons
    # exons_coords <- 
    #   getBM(attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand"),
    #         filters = c("ensembl_transcript_id"),
    #         values = transcript_id,
    #         mart = mart) %>% 
    #   dplyr::arrange(exon_chrom_start) %>%
    #   dplyr::mutate(dataset = str_replace(string = ensembl_dataset, pattern = "_gene_ensembl", replacement = ""), 
    #                 exon_index = ifelse(strand == "-1", n():1, 1:n())) 
    # 
    # # get coordinates and whole sequences of genes
    # genes_coords_seq <-
    #   getBM(attributes = c("start_position", "end_position", "gene_exon_intron"),
    #         filters = c("ensembl_transcript_id"),
    #         values = transcript_id,
    #         mart = mart) %>%
    #   dplyr::mutate(dataset = str_replace(string = ensembl_dataset, pattern = "_gene_ensembl", replacement = ""))
    # 
    # # get coordinates of introns
    # introns_coord <-
    #   exons_coords %>%
    #   dplyr::arrange(exon_chrom_start) %>%
    #   dplyr::mutate(int_chrom_start = exon_chrom_end + 1,
    #                 int_chrom_end = lead(exon_chrom_start) - 1,
    #                 exon_index = ifelse(strand == "-1", n():1, 1:n())) %>%
    #   dplyr::filter(!is.na(int_chrom_end)) %>%
    #   dplyr::mutate(intron_index = ifelse(strand == -1, exon_index - 1, exon_index + 1),
    #                 intron_index = ifelse(strand == -1, str_c(intron_index, "_", exon_index), str_c(exon_index, "_", intron_index))) %>%
    #   dplyr::select(-exon_chrom_start, -exon_chrom_end, -exon_index) %>%
    #   dplyr::select(chromosome_name, int_chrom_start, int_chrom_end, everything()) %>%
    #   dplyr::left_join(., tibble(dataset = animals,
    #                              animal_name = names(animals)),
    #                    by = "dataset") 
    # 
    # return(list(exons_coords, introns_coord, coding_seq, genes_coords_seq))
    # 
  })

coding_seq <- 
  lapply(bm_data, function(x) x[[3]]) %>% 
  dplyr::bind_rows(.) %>% 
  as.tibble(.)

coding_seq %>% 
  dplyr::filter(dataset == "mmusculus") %$%
  cdna %>% 
  nchar

exon_df <-
  lapply(bm_data, function(x) x[[1]]) %>%
  lapply(., function(x) mutate_all(x, as.character)) %>%
  dplyr::bind_rows(.) %>% 
  dplyr::left_join(., tibble(dataset = animals,
                             animal_name = names(animals)),
                   by = "dataset")  %>% 
  dplyr::mutate(liftOver = str_c(str_c("chr", chromosome_name), exon_chrom_start, exon_chrom_end, str_c("exon_", exon_index), sep = " "))

exon_df %>% 
  dplyr::filter(dataset == "mmusculus") %>%
  dplyr::select(liftOver) %>% 
  print(., row.names = FALSE)

introns_df <-
  lapply(bm_data, function(x) x[[2]]) %>%
  lapply(., function(x) mutate_all(x, as.character)) %>%
  dplyr::bind_rows(.) %>% 
  dplyr::filter(intron_index == "1_2")

# genes_coords_df <- 
#   lapply(bm_data, function(x) x[[2]]) %>% 
#   dplyr::bind_rows(.) %>% 
#   dplyr::as_tibble(.) %>% 
#   dplyr::mutate(n_char = nchar(gene_exon_intron)) %>% 
#   dplyr::filter(dataset == "ocuniculus") %$% 
#   gene_exon_intron

biomaRt::getSequence

### get sequences of chosen intron in gene from all mouse strains
# set intron in form LeftExonNumber_RightExonNumber (in sense direction)
intron_filt <- "1_2"

# get coordinates of introns from coordinates of exons
introns_coord <- 
  lapply(bm_data, function(x) x[[1]]) %>% 
  dplyr::bind_rows(.) %>% 
  group_by(dataset, ensembl_transcript_id) %>% 
  dplyr::select(chromosome_name, exon_chrom_start, exon_chrom_end, strand, dataset, ensembl_transcript_id) %>%
  dplyr::arrange(ensembl_transcript_id, exon_chrom_start) %>% 
  dplyr::mutate(start_int = exon_chrom_end + 1,
                end_int = lead(exon_chrom_start) - 1, 
                exon_index = ifelse(strand == "-1", n():1, 1:n())) %>%  
  dplyr::filter(!is.na(end_int)) %>%
  dplyr::mutate(intron_index = ifelse(strand == -1, exon_index - 1, exon_index + 1),
                intron_index = ifelse(strand == -1, str_c(intron_index, "_", exon_index), str_c(exon_index, "_", intron_index))) %>%
  dplyr::ungroup() %>% 
  dplyr::select(-exon_chrom_start, -exon_chrom_end, -exon_index, -ensembl_transcript_id) %>% 
  dplyr::select(chromosome_name, start = start_int, end = end_int, everything()) %>% 
  dplyr::distinct(.) %>% 
  dplyr::left_join(., tibble(dataset = animals,
                             animal_name = names(animals)), 
                   by = "dataset") %>% 
  dplyr::filter(intron_index == intron_filt)

# get coordinates and sequences of genes
genes_coords_seq <- 
  lapply(bm_data, function(x) x[[2]]) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::select(gene_start = start_position, gene_end = end_position, gene_exon_intron, dataset)

# join with intron coordinates with gene coordinates, get coordinates, subset gene sequences to chosen intron sequences
seq_data <- 
  dplyr::left_join(introns_coord, genes_coords_seq, by = "dataset") %>% 
  dplyr::mutate(subset_start = ifelse(strand == "-1", gene_end - end + 1, start - gene_start + 2), 
                subset_end = subset_start + (end - start - 1), 
                gene_exon_intron = stringr::str_sub(string = gene_exon_intron, start = subset_start, end = subset_end)) %>% 
  dplyr::select(-c(gene_start, gene_end, subset_start, subset_end)) 


### write sequences
# create dirs
outdir <- file.path(outpath, gene_symbol)
dir.create(outdir, showWarnings = F)

# get sequences
ensembl_seq <- 
  seq_data %>% 
  dplyr::filter(gene_exon_intron != "Sequence unavailable") %$% 
  gene_exon_intron %>% 
  DNAStringSet(.) 

# set names, order
names(ensembl_seq) <- 
  str_replace(seq_data$dataset, "(_|\\.).+", "") %>% 
  str_c(., gene_symbol, str_c("intron", seq_data$intron_index), sep = '_')

# write sequences
ensembl_seq %<>% 
  .[order(width(.))] %>% 
  .[(letterFrequency(., letters = "N", as.prob = T) %>% as.vector(.)) < 1] %T>%
  writeXStringSet(., file.path(outdir, str_c(gene_symbol, "mouse_strains", str_c("intron", intron_filt), length(.), "fa", sep = ".")), format = "fasta")

