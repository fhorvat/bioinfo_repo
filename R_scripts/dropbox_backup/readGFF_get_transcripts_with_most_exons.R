#!/home/students/fhorvat/R/bin/Rscript
### INFO: read GFF file, for each gene take transcipt with most exons
### DATE: 22. 5. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY

######################################################## LIBRARIES
library(stringr)
library(tibble)
library(tidyr)
library(data.table)

library(readr)
library(magrittr)
library(dplyr)
library(stringr)
library(tibble)
library(GenomicRanges)

######################################################## PATH VARIABLES
gtf_path_ENSEMBL <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Ensembl_GRCm38.86.20161128.gtf.gz"

######################################################## SOURCE FILES

######################################################## FUNCTIONS
GffToGRanges <- function(gff, filter = NULL){
  
  gff %<>%
    set_colnames(., c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")) %>%
    mutate(seqnames = ifelse(!str_detect(seqnames, "chr"), str_c("chr", seqnames), seqnames),
           seqnames = replace(seqnames, seqnames == "chrMT", "chrM")) %>%
    dplyr::filter(!str_detect(seqnames, "NT|GL|JH"))
  
  # gff has to have 9 columns
  if(ncol(gff) != 9)
    stop("Number of columns does not match gff format")
  
  # width of elements have to be positive
  if(any(as.integer(gff$end) < as.integer(gff$start))){
    warning("gff file contains ranges with negative widths...")
    gff <- gff[as.integer(gff$end) > as.integer(gff$start), ]
  }
  
  # filter - for example: exon/transcript/gene/CDS...
  if(!is.null(filter)){
    if(filter %in% gff$feature){
      gff = gff[gff$feature == filter, ]
    }else{
      stop("The given feature is not present in the gff file")
    }
  }
  
  ### get exon/transcript/gene ID for the 9th column of gff table
  # split stings in 9th column, split list again by length (not all genes have all info in that column)
  s <- str_split(string = gff$attribute, pattern = ";")
  z <- sapply(s, length)
  a <- split(s, z)
  gff <- gff[order(z), ]
  
  # create tables from 9th column with different number of columns in a list
  l <- lapply(a, function(x){
    
    d <- str_trim(str_replace(unlist(x, use.names = F), "^.+? ", "")) %>%
      str_replace_all(., "\"", "")
    
    m <- matrix(d, ncol = length(x[[1]]), byrow = T) %>%
      set_colnames(str_replace(str_trim(x[[1]]), " .+$", "")) %>%
      as_tibble()
    
    return(m)
    
  })
  
  # bind all tables from list in one with filling empty data with NA
  ids <-
    dplyr::bind_rows(l) %>%
    dplyr::select(matches("gene_id|transcript_id|exon_id|gene_biotype"))
  
  # set all strands which are not +/- to *
  gff$strand[!gff$strand %in% c('+', '-')] <- "*"
  
  # create and output GRanges
  granges <- GRanges(seqnames = gff$seqnames,
                     IRanges(as.integer(gff$start), as.integer(gff$end)),
                     strand = gff$strand,
                     frame = gff$frame,
                     feature.type = gff$feature,
                     .id = 1:nrow(gff))
  
  values(granges) <- cbind(values(granges), DataFrame(ids)[granges$.id, ])
  values(granges)$.id <- NULL
  
  return(granges)
}

######################################################## READ DATA
# read .gtf
gtf <-
  read_delim(file = gtf_path_ENSEMBL, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c")) %>% 
  GffToGRanges(., filter = "exon") 

######################################################## MAIN CODE
# create ENSEMBL transcripts ranges, takes one transcript for each gene (the one with the most exons)
gids <- unique(values(gtf)[c("gene_id", "transcript_id")])
gtf_trans <- split(gtf, gtf$transcript_id)
gtf_trans <- gtf_trans[order(-elementNROWS(gtf_trans))]
gtf_trans <- gtf_trans[!duplicated(gids$gene_id[match(names(gtf_trans), gids$transcript_id)])]

# take transcripts with single exon
gtf_trans_single <- 
  gtf_trans[elementNROWS(gtf_trans) == 1] %>% 
  unlist(.)
gtf_trans_single$ex.num <- 1
gtf_trans_single$ex.tot <- 1
gtf_trans_single <- split(gtf_trans_single, gtf_trans_single$transcript_id)

# take transcripts with more than one exon, count and enumerate exons in each transcript
gtf_trans <- gtf_trans[elementNROWS(gtf_trans) > 1]
gtf_trans <- unlist(gtf_trans)
d_val <- data.table(as.data.frame(values(gtf_trans)))
d_val$strand <- as.character(strand(gtf_trans))
d_val[d_val$strand == "+" , `:=`( COUNT = .N , IDX = 1:.N ) , by = transcript_id[strand == "+"]]
d_val[d_val$strand == "-" , `:=`( COUNT = .N , IDX = .N:.1) , by = transcript_id[strand == "-"]]
gtf_trans$ex.num <- d_val$IDX
gtf_trans$ex.tot <- d_val$COUNT
gtf_trans <- split(gtf_trans, as.character(gtf_trans$transcript_id))

# unite transcripts with one and more than one exons
gtf_trans <- 
  c(gtf_trans, gtf_trans_single) %>% 
  unlist() %>% 
  set_names(., NULL)
