library(dplyr)
library(stringr)
library(magrittr)
library(tibble)
library(GenomicRanges)
library(readr)

ensembl_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.gtf.gz"
ensembl_gtf <- readr::read_delim(file = ensembl_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c")) 
gff <- ensembl_gtf

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
  s <- stringr::str_split(string = gff$attribute, pattern = "; ")
  z <- sapply(s, length)
  a <- split(s, z)
  gff <- gff[order(z), ]
  
  # create tables from 9th column with different number of columns in a list
  l <- lapply(a, function(x){
    
    m <- 
      stringr::str_replace_all(unlist(x, use.names = F), "^.+? |\"|;", "") %>%
      matrix(., ncol = length(x[[1]]), byrow = T) %>%
      magrittr::set_colnames(str_replace(x[[1]], " .+$", "")) %>%
      tibble::as_tibble(.)
    
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
  
}
