### INFO: counts small RNA reads mapped to different categories in repeatMasker and in ENSEMBL annotation (hierarchically)
### DATE: Fri May 25 08:37:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# mouse
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Analysis/class_reads/ES_DcrTrans_2012.SOLiD")
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS
# class alignments hierarchically 
classReadsFactor <- function(bam_path, subject_ranges, yield = 1000000, isFirstInPair = NA){
  
  # get number of alignments in bam file
  read_number <-
    Rsamtools::countBam(file = bam_path, param = ScanBamParam(flag = scanBamFlag(isFirstMateRead = isFirstInPair))) %$%
    records
  
  # initialize alignment class vector
  align_class <- rep("placeholder", read_number)
  names(align_class) <- rep("nameplaceholder", read_number)
  
  # open connection to bam file in chunks
  bamfile <- BamFile(bam_path, yieldSize = yield)
  open(bamfile)
  
  # load chunks of alignments from bam file and classify each alignment
  while(length(chunk <- readGAlignmentsList(bamfile, param = ScanBamParam(what = "qname", flag = scanBamFlag(isFirstMateRead = isFirstInPair))))) {
    
    # unlist, set names of reads, transform to grglist (which gets ranges of only alignment part)
    chunk <- unlist(chunk)
    names(chunk) <- mcols(chunk)$qname
    chunk <- 
      GenomicRanges::grglist(chunk) %>% 
      unlist(.)
    
    # find overlaps between two GRanges - chunk of alignments and annotation
    hits <- findOverlaps(chunk, subject_ranges, ignore.strand = T)
    
    # class if there are any hits
    if(length(hits) > 0){
      
      # get hits in alignments 
      read_hits <- 
        extractList(chunk, as(queryHits(hits), "List")) %>% 
        unlist(.)
      
      # get hits in annotation
      subject_hits <-
        extractList(count_ranges, as(subjectHits(hits), "List")) %>% 
        unlist(.)
      
      # create vector of gene biotypes
      sense_strand <- ifelse(strand(read_hits) == strand(subject_hits), "sense", "antisense")
      gene_biotype_vector <- subject_hits$gene_biotype
      gene_biotype_vector[gene_biotype_vector == "miRNA"] <- str_c(gene_biotype_vector[gene_biotype_vector == "miRNA"], ".", sense_strand[gene_biotype_vector == "miRNA"])
      gene_biotype_vector[gene_biotype_vector == "protein_coding"] <- str_c(gene_biotype_vector[gene_biotype_vector == "protein_coding"], ".", sense_strand[gene_biotype_vector == "protein_coding"])
      gene_biotype_vector <- replace(gene_biotype_vector, 
                                     !(gene_biotype_vector %in% class_hier[class_hier != "other"]), 
                                     "other")
      
      # class each unique read ID to category
      read_class_vector <-
        factor(x = gene_biotype_vector, levels = class_hier) %>%
        as.numeric(.) %>%
        split(., factor(names(read_hits))) %>%
        purrr::map(., min) %>%
        unlist(.)
      
      # find first element which is not filled
      last_element <-
        which(align_class == "placeholder") %>%
        min(.)
      
      # set range in intialized vector
      class_range <- c(last_element, (last_element - 1 + length(read_class_vector)))
      
      # create named vector (read names and class as name)
      align_class[class_range[1]:class_range[2]] <- names(read_class_vector)
      names(align_class)[class_range[1]:class_range[2]] <- class_hier[read_class_vector]
      
      # filter reads out
      chunk <- chunk[!(names(chunk) %in% names(read_hits))]
      
    }
    
    # set alignements which are not overlaping any category as "not_annotated"
    if(length(chunk) > 0){
      
      # find first element which is not filled
      last_element <-
        which(align_class == "placeholder") %>%
        min(.)
      
      # set range in intialized vector
      class_range <- c(last_element, (last_element - 1 + length(unique(names(chunk)))))
      
      # add to named vector (read names and class as name)
      align_class[class_range[1]:class_range[2]] <- unique(names(chunk))
      names(align_class)[class_range[1]:class_range[2]] <- "not_annotated"
      
    }
    
  }
  
  # close connection to .bam
  close(bamfile)
  
  # remove placeholders from vector
  align_class <- align_class[align_class != "placeholder"]
  
  # return vector with read names
  return(align_class)
  
}

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# get gene info path
info_path <- list.files(genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.geneInfo.csv", full.names = T)

# get paths of reduced exons
exons_path <- list.files(genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.reducedExons.RDS", full.names = T)

# get repeatMasker path
rmsk_path <- list.files(genome_path, pattern = "rmsk.*[0-9]{6}.clean.fa.out.gz", full.names = T)

# get bam file path, name, experiment
bam_path <- "%BAM_PATH"
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Mapped/STAR_mm10.pCAG_EGFP_MosIR/s_RS7_NC_r1.SE.genome.Aligned.sortedByCoord.out.bam"
bam_name <- basename(bam_path) %>% str_remove_all(., ".SE.genome.Aligned.sortedByCoord.out.bam|.SE|.bam")
experiment_name <- "ES_DcrTrans_2012"

######################################################## READ DATA
# read info about genes
genes_info <- readr::read_csv(info_path)

# gtf exons by genes
exons_gr <- readRDS(file = exons_path)

# read repeatMasker
rmsk_df <- readr::read_delim(rmsk_path, delim = "\t")

######################################################## MAIN CODE
### set class hierarchy
class_hier <- c("Mos_mRNA", "miRNA.sense", "protein_coding.sense", "rRNA", "SINE", "LINE", "LTR", "other_repeat", "annotated_pseudogene", "other", "not_annotated")

### prepare ranges to count over - repeatMasker and exons of genes from ENSEMBL annotation
# get GRanges from repeatMasker
rmsk_gr <- 
  rmsk_df %>% 
  dplyr::mutate(gene_id = str_c(seqnames, ":", start, "-", end, "|", strand, "|", repName)) %>% 
  dplyr::select(seqnames:strand, gene_id, gene_biotype = repClass) %>% 
  dplyr::mutate(gene_biotype = replace(x = gene_biotype, 
                                       list = !(gene_biotype %in% c("rRNA", "LINE", "SINE", "LTR")), 
                                       values = "other_repeat")) %>% 
  GenomicRanges::GRanges(.)

# prepare exons
exons_gr_clean <-
  exons_gr %>% 
  as.data.frame(.) %>% 
  tibble::as.tibble(.) %>% 
  dplyr::left_join(., genes_info %>% dplyr::select(gene_id, gene_biotype), by = "gene_id") %>% 
  dplyr::select(seqnames:strand, gene_id, gene_biotype) %>% 
  dplyr::mutate(gene_biotype = 
                  replace(gene_biotype, gene_biotype == "Mt_rRNA", "rRNA") %>% 
                  replace(., str_detect(., "pseudogene"), "annotated_pseudogene")) %>%
  GenomicRanges::GRanges(.)

# mos plasmid
mos_gr <- GenomicRanges::GRanges(seqnames = "pCAG-EGFP_MosIR", 
                                 ranges = IRanges(start = 1, end = 6660), 
                                 gene_id = "pCAG-EGFP_MosIR", 
                                 gene_biotype = "Mos_mRNA")

# join repeatMasker and ENSEMBL
count_ranges <- c(rmsk_gr, exons_gr_clean, mos_gr)

### class reads and save table
reads_class <- classReadsFactor(bam_path = bam_path, 
                                subject_ranges = count_ranges, 
                                yield = 1000000, 
                                isFirstInPair = NA)

# summarize read classification
read_class_vector <-
  factor(x = names(reads_class), levels = class_hier) %>%
  as.numeric(.) %>%
  split(., factor(reads_class)) %>%
  purrr::map(., min) %>%
  unlist(.)

# sum sense and antisense read groups
read_hits_final <- 
  tibble::tibble(read_id = names(read_class_vector), 
                 read_group = class_hier[read_class_vector]) %>% 
  unique(.) %>% 
  dplyr::group_by(read_group) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::mutate(read_group = class_hier[read_group], 
                sample = bam_name, 
                experiment = experiment_name) %T>% 
  readr::write_csv(x = ., path = file.path(outpath, str_c("smallRNA", experiment_name, bam_name, "csv", sep = ".")))
