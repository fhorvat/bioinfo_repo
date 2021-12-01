### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd(".")

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

library(biomaRt)
library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/human/hg38.GRCh38.GCA_000001405.15"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames\\.geneInfo\\.csv$"), full.names = T)

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

######################################################## MAIN CODE
# Aravin et al. (10.1371/journal.pgen.1000764): 
# MILI-TDRD1 - The first type of granules, pi-bodies, contains the MILI-TDRD1 module of the piRNA pathway
# MIWI2-TDRD9-MAEL - The second type of granules, piP-bodies, harbors the MIWI2-TDRD9-MAEL module of the piRNA 
# MVH - mouse VASA homolog (MVH) protein, an RNA helicase
# 
# Ernst et al. (10.1038/s41467-017-01049-7): 
# Maelstrom (MAEL) - piRNA precursor molecules are then exported into the cytoplasm, which might be facilitated by Maelstrom (MAEL)
# PLD6 (MITOPLD) - the endonuclease that generates the 5' end of primary piRNAs = Zuccchini
# Tdrkh - The process is facilitated by Tdrkh, a Tudor and KH domain-containing protein, whose absence causes the accumulation of 31–37 nucleotide (nt) long intermediates
# HENMT1 - Mature 3' ends are then 2'-O methylated by the RNA methyltransferase HENMT1
# RNF17 - The disengagement of MIWI from the ping-pong cycle was shown to be due to active repression of the ping-pong cycle by RNF17, a Tudor family member, in meiotic cells
# A-MYB - Transcription of these clusters is facilitated by the ancestral transcription factor A-MYB
# DDX4 - Vasa (also known as Ddx4) is an ATP-dependent RNA helicase

# # Pandey et al. (10.1073/pnas.1316316110)
# TDRD12 - Tudor domain containing 12 (TDRD12) is essential for secondary PIWI interacting RNA biogenesis in mice
# 
# Vagin et al. (10.1101/gad.1814809): 
# Piwi proteins were found in complex with PRMT5/WDR77, an enzyme that dimethylates arginine residues
# 
# Chen et al. (10.1073/pnas.0911640106): 
# Tdrd1, Tdrd4/RNF17, and Tdrd6, Tdrkh/Tdrd2
# In this analysis, we observed several previously known Piwi-interacting proteins or piRNA pathway components such as Ddx4, Kif17, Fxr1, and Mov10l1

# Siomi et al. (10.1101/gad.1899210)
# PIWI proteins have been shown recently to contain symmetrical dimethyl arginines (sDMAs), 
# and this modification is mediated by the methyltransferase PRMT5 (also known as Dart5 or Capsuleen)

### piRNA-pathway associated proteins
# create table
pirna_genes_tb <- 
  tibble(gene_name = c("Piwil2", "Piwil1", "Piwil4", "Piwil3",
                       "Mov10l1", "Ddx4", "Mael", 
                       "Pld6", "Henmt1",
                       "Tdrd1", "Tdrkh", "Rnf17", "Tdrd6", "Tdrd9", "Tdrd12", 
                       "Prmt5", "Wdr77", "Kif17"),
         other_name = c("MILI", "MIWI", "MIWI2", "MIWI3",
                        "MOV10L1", "DDX4/VASE", "MAEL/Maelstorm", 
                        "PLD6/MITOPOLD/Zucchini", "HENMT1", 
                        "TDRD1", "TDRD2/TDRKH", "TDRD4/RNF17", "TDRD6", "TDRD9", "TDRD12", 
                        "PRMT5/Dart5/Capsuleen", "WDR77", "KIF17"))

# filter genes info to get gene IDs
gene_names_tb <- 
  genes_info %>% 
  dplyr::filter(gene_name %in% toupper(pirna_genes_tb$gene_name)) %>% 
  tidyr::unite(col = "coordinates", seqnames, start, end, sep = " ")


### download gene info from Ensembl
# Ensembl versions
ensembl_url <-
  tibble(ens_version = c(100, 99, 98, 96, 95, 94, 93, 92, 91, 89, 86),
         date = c("Apr 2020", "Jan 2020", "Sep 2019", "Apr 2019", "Jan 2019", "Oct 2018", "Jul 2018", "Apr 2018", "Dec 2017", "May 2017", "Oct 2016"),
         URL_archive = c("www.ensembl.org",
                         "http://jan2020.archive.ensembl.org",
                         "http://sep2019.archive.ensembl.org",
                         "http://apr2019.archive.ensembl.org",
                         "http://jan2019.archive.ensembl.org",
                         "http://oct2018.archive.ensembl.org",
                         "http://jul2018.archive.ensembl.org",
                         "http://apr2018.archive.ensembl.org",
                         "http://dec2017.archive.ensembl.org",
                         "http://may2017.archive.ensembl.org",
                         "http://oct2016.archive.ensembl.org")) %>%
  dplyr::filter(ens_version == ensembl_version) %$%
  URL_archive

# load ENSEMBL mart
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = stringr::str_c("hsapiens", "_gene_ensembl"), host = ensembl_url)

# get gene IDs, transcript ID's and transcript names
gene_transcript_tb <- 
  getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "external_transcript_name"), 
        filters = "ensembl_gene_id",
        values = gene_names_tb$gene_id, 
        mart = mart) %>%
  tibble::as_tibble(.) %>% 
  dplyr::mutate(transcript_index = 
                  str_extract(external_transcript_name, "(?<=-)[0-9]+") %>% 
                  as.numeric(.)) %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::filter(transcript_index == min(transcript_index)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., gene_names_tb %>% dplyr::select(gene_id, gene_name), by = c("ensembl_gene_id" = "gene_id")) %>% 
  dplyr::select(gene_name, ensembl_transcript_id, ensembl_gene_id)

# get protein sequences
protein_tb <- 
  getSequence(id = gene_transcript_tb$ensembl_transcript_id, 
              type = "ensembl_transcript_id", 
              seqType = "peptide", 
              mart = mart) %>% 
  as_tibble(.) %>% 
  dplyr::left_join(., gene_transcript_tb, by = "ensembl_transcript_id") %>% 
  dplyr::select(gene_name, gene_id = ensembl_gene_id, transcript_id = ensembl_transcript_id, protein_seq = peptide) 


### save files
# create AAStringSet
protein_seq <- Biostrings::AAStringSet(protein_tb$protein_seq)
names(protein_seq) <- str_c(protein_tb$gene_name)

# save
writeXStringSet(protein_seq, filepath = file.path(outpath, str_c("ensembl.99.GRCh38.p12.20200505.piRNA_pathway_genes.AA.fasta")))

# save table
protein_tb_tidy <- 
  protein_tb %>% 
  dplyr::select(gene_id, transcript_id) %>% 
  dplyr::left_join(., gene_names_tb, by = "gene_id") %>% 
  dplyr::select(gene_name, gene_id:strand, gene_description)

# join with original table
pirna_genes_tb %>% 
  dplyr::mutate(gene_name = toupper(gene_name)) %>% 
  dplyr::left_join(., protein_tb_tidy, by = "gene_name") %T>% 
  readr::write_csv(., path = file.path(outpath, str_c("ensembl.99.GRCh38.p12.20200505.piRNA_pathway_genes.infoTable.csv")))
  
  
