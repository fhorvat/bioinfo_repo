### INFO: counts small RNA reads mapped to different categories in repeatMasker and in ENSEMBL annotation (hierarchically)
### DATE: Fri May 25 08:37:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
# wideScreen()

######################################################## WORKING DIRECTORY
# mouse
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

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)

######################################################## SOURCE FILES

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
  
  # return vector with read names
  return(align_class)
  
}

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# get gene info path
info_path <- list.files(genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.geneInfo.csv", full.names = T)

# get paths of reduced exons
exons_path <- list.files(genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.reducedExons.RDS", full.names = T)

# get repeatMasker path
rmsk_path <- list.files(genome_path, pattern = "rmsk.*[0-9]{6}.clean.fa.out.gz", full.names = T)

# get miRBase gff path
mirbase_path <- list.files(genome_path, pattern = "miRBase.*gff3", full.names = T)

# get path of inverted repeat coordinates from mRNA
ir_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/mm10.all_plasmids/IR_mRNAs.plasmid_coordinates/IR_mRNA.plasmid_coordinates.clean.csv"

# get bam file path, name, experiment and IR name
# bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/NIH3T3_transfected.2018/Data/Mapped/STAR_mm10.plasmids/s_pCag_EGFP_mosIR_r1.SE.mis_0.bam"
bam_path <- "%BAM_PATH"
bam_name <- basename(bam_path) %>% str_remove_all(., ".bam")
experiment_name <- str_remove(bam_path, "/Data.*") %>% basename(.)
ir_name <- str_remove_all(bam_name, "^s_|_r[1,2].SE.mis_[0-2]") %>% str_replace(., "pCag_EGFP", "pCag-EGFP")

######################################################## READ DATA
# read info about genes
genes_info <- readr::read_csv(info_path)

# gtf exons by genes
exons_gr <- readRDS(file = exons_path)

# read repeatMasker
rmsk_gr <- readr::read_delim(rmsk_path, delim = "\t")

# read miRBase gff
mirna_gr <- rtracklayer::import.gff(con = mirbase_path)

# read coordinates of inverted repeats from mRNA
ir_gr <- readr::read_csv(file = ir_path)

######################################################## MAIN CODE
### set class hierarchy
class_hier <- c(ir_name,  
                "miRNA.mature.sense", "miRNA.other.sense", "protein_coding.sense", 
                "rRNA", "SINE", "LINE", "LTR", "other_repeat", 
                "annotated_pseudogene", "other", "not_annotated")

### prepare ranges to count over - repeatMasker, exons of genes from ENSEMBL annotation, miRBase, mosIR plasmid mRNA
# get GRanges from repeatMasker
rmsk_gr %<>%
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

# prepare mRNA IR coordinates on plasmids
ir_gr %<>% GenomicRanges::GRanges(.)
mcols(ir_gr) <- mcols(ir_gr)[, c("arm")]
names(mcols(ir_gr)) <- "gene_id"
mcols(ir_gr)$gene_biotype <- as.character(seqnames(ir_gr))

# join repeatMasker, Ensembl, miRBase and IR GRanges
count_ranges <- c(rmsk_gr, exons_gr, mirna_gr, ir_gr)

### class reads and save table
reads_class <- classReadsFactor(bam_path = bam_path,
                                subject_ranges = count_ranges,
                                yield = 1000000,
                                isFirstInPair = NA)


### loop through unique names and sum read categories
# get unique names
reads_class_unique <- unique(reads_class)

# create sum table
reads_class_sum <- tibble(read_group = class_hier,
                          count = 0)

# set loop step
loop_step <- 100000

# set loop index
loop_index <- seq(0, length(reads_class_unique), by = loop_step)
if(loop_index[length(loop_index)] < length(reads_class_unique)){
  loop_index <- c(loop_index, length(reads_class_unique))
}

# loop
for(n in 1:(length(loop_index) - 1)){
  
  # subset vector
  reads_class_subset <- reads_class[reads_class %in% reads_class_unique[(loop_index[n] + 1):loop_index[n + 1]]]
  
  # summarize read classification
  read_class_vector <-
    factor(x = names(reads_class_subset), levels = class_hier) %>%
    as.numeric(.) %>%
    split(., factor(reads_class_subset)) %>%
    purrr::map(., min) %>%
    unlist(.)
  
  # summarize
  reads_class_sum_subset <-
    tibble::tibble(read_id = names(read_class_vector),
                   read_group = class_hier[read_class_vector]) %>%
    unique(.) %>%
    dplyr::group_by(read_group) %>%
    dplyr::summarise(count = n())
  
  # join with sum vector
  reads_class_sum %<>%
    dplyr::left_join(., reads_class_sum_subset, by = "read_group") %>%
    dplyr::mutate(count.y = replace(count.y, is.na(count.y), 0),
                  count = count.x + count.y) %>%
    dplyr::select(read_group, count)
  
}

# add info about sample and experiment, save
reads_class_sum %<>%
  dplyr::mutate(sample_id = bam_name,
                experiment = experiment_name) %T>%
  readr::write_delim(x = ., path = file.path(outpath, str_c("read_class", bam_name, "txt", sep = ".")), delim = "\t")