### INFO: counts small RNA reads mapped to different categories in repeatMasker and in ENSEMBL annotation (hierarchically)
### DATE: Fri May 25 08:37:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# mouse
setwd(".")
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
# lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
# source(file.path(lib_path, "wideScreen.R"))
# source(file.path(lib_path, "headt.R"))
# source(file.path(lib_path, "asdf.R"))
# wideScreen()

######################################################## FUNCTIONS
# class alignments hierarchically 
classReadsFactor <- function(bam_path, subject_ranges, yield = 1000000, isFirstInPair = NA, nm_tag, nm_tag_value){
  
  # set tag filter
  tag_filter <- list(nm_tag_value) %>% magrittr::set_names(nm_tag)
  
  # get number of alignments in bam file
  read_number <-
    Rsamtools::countBam(file = bam_path, 
                        param = ScanBamParam(flag = scanBamFlag(isFirstMateRead = isFirstInPair), 
                                             tagFilter = tag_filter)) %$%
    records
  
  # initialize alignment class vector
  align_class <- rep("placeholder", read_number)
  names(align_class) <- rep("nameplaceholder", read_number)
  
  # open connection to bam file in chunks
  bamfile <- BamFile(bam_path, yieldSize = yield)
  open(bamfile)
  
  # load chunks of alignments from bam file and classify each alignment
  while(length(chunk <- readGAlignmentsList(file = bamfile, 
                                            param = ScanBamParam(what = "qname", 
                                                                 flag = scanBamFlag(isFirstMateRead = isFirstInPair), 
                                                                 tagFilter = tag_filter)))) {
    
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
      gene_biotype_vector[str_detect(gene_biotype_vector, "miRNA")] <- str_c(gene_biotype_vector[str_detect(gene_biotype_vector, "miRNA")], ".", 
                                                                             sense_strand[str_detect(gene_biotype_vector, "miRNA")])
      gene_biotype_vector[gene_biotype_vector == "protein_coding"] <- str_c(gene_biotype_vector[gene_biotype_vector == "protein_coding"], ".", 
                                                                            sense_strand[gene_biotype_vector == "protein_coding"])
      gene_biotype_vector <- replace(gene_biotype_vector, 
                                     !(gene_biotype_vector %in% class_hier[class_hier != "other"]), 
                                     "other")
      
      # class each unique read ID to category
      read_class_vector <-
        factor(x = gene_biotype_vector,
               levels = class_hier) %>%
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
  
  
  ### loop through unique names and sum read categories
  # get unique names
  align_class_unique <- unique(align_class)
  
  # create sum table
  align_class_sum <- tibble(read_group = class_hier, count = 0) 
  
  # set loop step
  loop_step <- yield
  
  # set loop index
  loop_index <- seq(0, length(align_class_unique), by = loop_step)
  if(loop_index[length(loop_index)] < length(align_class_unique)){
    loop_index <- c(loop_index, length(align_class_unique))
  }
  
  # loop
  for(n in 1:(length(loop_index) - 1)){
    
    # subset vector
    align_class_subset <- align_class[align_class %in% align_class_unique[(loop_index[n] + 1):loop_index[n + 1]]]
    
    # summarize read classification
    read_class_vector <-
      factor(x = names(align_class_subset), levels = class_hier) %>%
      as.numeric(.) %>%
      split(., factor(align_class_subset)) %>%
      purrr::map(., min) %>%
      unlist(.)
    
    # summarize
    align_class_sum_subset <-
      tibble::tibble(read_id = names(read_class_vector), 
                     read_group = class_hier[read_class_vector]) %>% 
      unique(.) %>%
      dplyr::group_by(read_group) %>%
      dplyr::summarise(count = n())
    
    # join with sum vector
    align_class_sum %<>%  
      dplyr::left_join(., align_class_sum_subset, by = "read_group") %>% 
      dplyr::mutate(count.y = replace(count.y, is.na(count.y), 0), 
                    count = count.x + count.y) %>% 
      dplyr::select(read_group, count)
    
  }
  
  # add mismatch value to tibble
  align_class_sum %<>%
    dplyr::mutate(mismatch = nm_tag_value)
  
  # return vector with read names
  return(align_class_sum)
  
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

# get miRBase gff path
mirbase_path <- list.files(genome_path, pattern = "miRBase.*gff3", full.names = T)

# get .BED with mosIR coordinates
mosir_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/pCAG-EGFP_MosIR.coordinates.bed"

# get bam file path, name, experiment
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/ES_DcrTrans_2012/Data/Mapped/Shrimp_mm10/s_ES_Dcr_Rsc_utran_r2.SE.bam"
# bam_path <- "%BAM_PATH"
bam_name <- basename(bam_path) %>% str_remove_all(., ".SE.genome.Aligned.sortedByCoord.out.bam|.SE|.bam")
experiment_name <- str_remove(bam_path, "/Data.*") %>% basename(.)

######################################################## READ DATA
# read info about genes
genes_info <- readr::read_csv(info_path)

# gtf exons by genes
exons_gr <- readRDS(file = exons_path)

# read repeatMasker
rmsk_df <- readr::read_delim(rmsk_path, delim = "\t")

# read miRBase gff
mirna_gr <- rtracklayer::import.gff(con = mirbase_path) 

# read coordinates of mosIR on plasmid
mosir_gr <- rtracklayer::import.bed(con = mosir_path)

######################################################## MAIN CODE
### set class hierarchy
class_hier <- c("Mos_mRNA", "miRNA.mature.sense", "miRNA.other.sense", "protein_coding.sense", "rRNA", "SINE", "LINE", "LTR", "other_repeat", "annotated_pseudogene", "other", "not_annotated")

### prepare ranges to count over - repeatMasker, exons of genes from ENSEMBL annotation, miRBase, mosIR plasmid mRNA
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
exons_gr %<>% 
  as.data.frame(.) %>% 
  tibble::as.tibble(.) %>% 
  dplyr::left_join(., genes_info %>% dplyr::select(gene_id, gene_biotype), by = "gene_id") %>% 
  dplyr::select(seqnames:strand, gene_id, gene_biotype) %>% 
  dplyr::mutate(gene_biotype = 
                  replace(gene_biotype, gene_biotype == "Mt_rRNA", "rRNA") %>% 
                  replace(., str_detect(., "pseudogene"), "annotated_pseudogene") %>% 
                  replace(., . == "miRNA", "miRNA.other")) %>%
  GenomicRanges::GRanges(.)

# prepare  ranges of mature miRNA
mirna_gr <- mirna_gr[mcols(mirna_gr)$type == "miRNA"]
mcols(mirna_gr) <- mcols(mirna_gr)[, c("Name")]
names(mcols(mirna_gr)) <- "gene_id"
mcols(mirna_gr)$gene_biotype <- "miRNA.mature"

# prepare MosIR coordinates 
names(mcols(mosir_gr)) <- "gene_id" 
mcols(mosir_gr)$gene_biotype <- "Mos_mRNA"

# join repeatMasker and ENSEMBL
count_ranges <- c(rmsk_gr, exons_gr, mirna_gr, mosir_gr)

# class reads with 0, 1 or 2 mismatches
reads_class <- purrr::map(.x = c(0, 1, 2), .f = function(mismatch){
  classReadsFactor(bam_path = bam_path, 
                   subject_ranges = count_ranges, 
                   yield = 1000000, 
                   isFirstInPair = NA, 
                   nm_tag = "NM", 
                   nm_tag_value = mismatch) 
})
                          

# add info about sample and experiment, save
reads_class_all <- 
  dplyr::bind_rows(reads_class) %>% 
  dplyr::mutate(sample = bam_name, 
                experiment = experiment_name) %T>% 
  readr::write_delim(x = ., path = file.path(outpath, str_c("read_class", bam_name, "txt", sep = ".")), delim = "\t")

